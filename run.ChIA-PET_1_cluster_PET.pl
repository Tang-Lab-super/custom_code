#!/usr/bin/perl
use strict;
use warnings;
use threads;
use threads::shared;
use Data::Dumper;
use FindBin qw($Bin);
#use lib "$Bin/../lib";
use lib "/data/public/software/CHIAPET_PIPELINE/lib/";
use GenomeSize;

## fix bug on May 9th, 2017
## $cpu_num didn't pass to subrotinate 'Run_Clustering_Multi_Threads'

my $input_bedpe          = shift; ## bedpe is indexed by pairix
my $output_name          = shift;
my $ref_genome           = shift;
my $self_ligation_cutoff = shift; ## span cutoff for self-ligation
my $large_span_cutoff    = shift; ## exclude large span PET for clustering; RNAPII data (2 Mb)
my $clustering_type      = shift; ## all, intra or inter
my $cpu_num              = shift;
my $CMD_pairix           = shift; # pairix (https://github.com/4dn-dcic/pairix); requery PET in a fast way
my $ext                  = shift; # extension from one side of a fragment

####################### default extension for clustering ###################################
#my $ext = 500;

unless ( -e $CMD_pairix ){ print STDERR 'Command doesn\'t exist:', "\n\t", $CMD_pairix, "\n"; exit 1; }


###################### reference genome chromosome index ##################################
my $genome    = GenomeSize->new( genome => $ref_genome );
my $chr_index = $genome->Get_chr_index();
my @chr_names = sort { $chr_index->{ $a } <=> $chr_index->{ $b } } keys %{ $chr_index };

my $prefix = Random_String();

########################### clustering with multi-threads ##################################
print STDERR "Cluster PETs using multi-threads\n";
my $output_clusters_files = Run_Clustering_Multi_Threads ( \@chr_names, $clustering_type, $input_bedpe, $prefix, $self_ligation_cutoff, $large_span_cutoff, $ext, $cpu_num, $CMD_pairix );

print STDERR "Merge outputs\n";
&Merge_Files ( $output_clusters_files, $output_name );





##################################################################################################################################################
############################### subrotinates ############################################

sub Run_Clustering_Multi_Threads {
    my ( $chr_names, $clustering_type, $indexed_bedpe, $prefix, $self_ligation_cutoff, $large_span_cutoff, $ext, $cpu_num, $CMD_pairix ) = @_;

    my $chr_pairs      = Create_Chr_Pairs ( $chr_names, $clustering_type );
    my @chr_pairs_list = @{ $chr_pairs };

    my @output_clusters :shared;
    my $output_clusters_ref = \@output_clusters;


## require all numbers of CPU
    if( scalar @chr_pairs_list >= $cpu_num ){
	while( (scalar @chr_pairs_list - $cpu_num) >= 0 ){
	    my @threads;
	    my @sub_group = splice( @chr_pairs_list, 0, $cpu_num );

	    foreach my $pair ( @sub_group ){
		my $thr = threads->create( 'Run_Clustering_PET', 
					   ( $pair->[0], $pair->[1], $input_bedpe, $prefix, $self_ligation_cutoff, $large_span_cutoff, $output_clusters_ref, $ext, $CMD_pairix ) );
		push @threads, $thr;
	    }
	    foreach ( @threads ){
		my $tid = $_->tid();
		print STDERR "Thread collected: $tid\n";

		$_->join();
		
		if ( @chr_pairs_list ){
		    my $pair_add = shift @chr_pairs_list;
		    my $thr_add = threads->create( 'Run_Clustering_PET',
						   ( $pair_add->[0], $pair_add->[1], $input_bedpe, $prefix, $self_ligation_cutoff, $large_span_cutoff, $output_clusters_ref, $ext, $CMD_pairix ) );
		    push @threads, $thr_add;
		}
	    }
	}
    }

## needs few numbers of CPU
    else {
	my @threads;
	for my $pair ( @chr_pairs_list ){
	    my $thr = threads->create( 'Run_Clustering_PET', 
				       ( $pair->[0], $pair->[1], $input_bedpe, $prefix, $self_ligation_cutoff, $large_span_cutoff, $output_clusters_ref, $ext, $CMD_pairix ) );
	    push @threads, $thr;
	}
	foreach ( @threads ){
	    my $tid = $_->tid();
	    print STDERR "Thread collated: $tid\n";

	    $_->join();
	}
    }
    
    return $output_clusters_ref;
}


sub Run_Clustering_PET {
    my ( $chr_1, $chr_2, $indexed_bedpe, $prefix, $self_ligation_cutoff, $large_span_cutoff, $output_clusters_ref, $ext, $CMD_pairix ) = @_;

    my $output_clusters = $prefix.'.'.$chr_1.'_'.$chr_2.'.clusters.bedpe';
    push @{ $output_clusters_ref }, $output_clusters;

    open( my $fh_out, ">", $output_clusters ) or die $!;

    my $pets = Parsing_Indexed_Bedpe ( $indexed_bedpe, $chr_1, $chr_2, $self_ligation_cutoff, $large_span_cutoff, $CMD_pairix );

    &Clustering_PETs ( $pets, $chr_1, $chr_2, $ext, $fh_out );
    
    close $fh_out;

    return 0;
}


sub Clustering_PETs {
    my $pets     = shift;
    my $chr_head = shift;
    my $chr_tail = shift;
    my $ext      = shift;
    my $fh_out   = shift;

    my ( $head_start, $head_end, $head_overlap_flag, @head_overlap );

    my $i;
    for my $start_1 ( sort { $a <=> $b } keys %{ $pets } ){
	for my $end_1 ( sort { $a <=> $b } keys %{ $pets->{ $start_1 } } ){
	    $i++;
	    print STDERR "Processing $chr_head, $chr_tail: [$i]\r" unless $i % 10000;
	    

	    if ( ! $head_start && ! $head_end ){
		$head_start        = $start_1;
		$head_end          = $end_1;
		$head_overlap_flag = 0;
		
	    } else {
		my ( $head_overlap_s, $head_overlap_e ) = Overlapping ( $head_start, $head_end, $start_1, $end_1,  $ext ); # input in sorted order

#case 1
		if ( ! $head_overlap_s && ! $head_overlap_e && ! $head_overlap_flag ){

		    my @single_head;
		    push @single_head, [$head_start, $head_end];

		    &Clustering_and_output ( $pets, \@single_head, $chr_head, $chr_tail, $ext, $fh_out );

		    $head_start        = $start_1;
		    $head_end          = $end_1;
		    $head_overlap_flag = 0;

		} 
#case 2
		elsif ( $head_overlap_s && $head_overlap_e && ! $head_overlap_flag ) {

		    push @head_overlap, [$head_start, $head_end];

		    $head_start        = $start_1;
		    $head_end          = $end_1;
		    $head_overlap_flag = 1;

		} 
#case 3
		elsif ( ! $head_overlap_s && ! $head_overlap_e && $head_overlap_flag ){

		    push @head_overlap, [$head_start, $head_end];
		    &Clustering_and_output ( $pets, \@head_overlap, $chr_head, $chr_tail, $ext, $fh_out );

		    undef @head_overlap;

		    $head_start        = $start_1;
		    $head_end          = $end_1;
		    $head_overlap_flag = 0;
		}
#case 4
		elsif ( $head_overlap_s && $head_overlap_e && $head_overlap_flag ){
		    

		    push @head_overlap, [$head_start, $head_end];

		    $head_start        = $start_1;
		    $head_end          = $end_1;
		    $head_overlap_flag = 1;

		}
		else {
		    print STDERR 'Some case is missing: ', join ("\t", ($head_start, $head_end, $start_1, $end_1)), "\n";
		}
	    }
	} # end_1
    } # start_1


#case 5
# the last PET 
    if ( $head_overlap_flag ){

	push @head_overlap, [$head_start, $head_end];


	&Clustering_and_output ( $pets, \@head_overlap, $chr_head, $chr_tail, $ext, $fh_out );

	undef @head_overlap;

    } 
#    elsif ( $head_start && $head_end ) {
#	next;
#	my @tail_p = @{ $pets->{ $head_start }{ $head_end } };
#	for my $tail ( @tail_p ){
#	    my $tail_start = $tail->[0];
#	    my $tail_end   = $tail->[1];
#	    print $fh_out join ( "\t", ( $chr_head, $head_start, $head_end, $chr_tail, $tail_start, $tail_end, 1 ) ), "\n";
#	}
#    }
    print STDERR "\n";

    return 0;
}



sub Clustering_and_output {
    my $pets         = shift;
    my $head_overlap = shift;
    my $chr_head     = shift;
    my $chr_tail     = shift;
    my $ext          = shift;
    my $fh_out       = shift;

# get head overlapped PETs
#    my ( @head_overlap_pets, %clusters );
    my @head_overlap_pets;

    for my $head ( sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{ $head_overlap }  ){
	my $h_s = $head->[0];
	my $h_e = $head->[1];
	
	my @tail_pos = @{ $pets->{ $h_s }{ $h_e } };

	for my $tail ( @tail_pos ){
	    my $t_s = $tail->[0];
	    my $t_e = $tail->[1];
	    push @head_overlap_pets, join("-", ( $h_s, $h_e, $t_s, $t_e, 1 ));
	}
    }

    my ( @clustered_pets, $overlap_flag );

  PET:
    while ( my $pet_1 = shift @head_overlap_pets ){

	next if ( $pet_1 eq 'NA' );

	my ( $h_s_1, $h_e_1, $t_s_1, $t_e_1, $count_1 ) = split(/-/, $pet_1);

	for ( my $i = 0; $i <= $#head_overlap_pets; $i++ ){

	    next if ( $head_overlap_pets[$i] eq 'NA' );

	    my $pet_2 = $head_overlap_pets[$i];
	    my ( $h_s_2, $h_e_2, $t_s_2, $t_e_2, $count_2 ) = split(/-/, $pet_2);

	    my ( $head_o_s, $head_o_e ) = Overlapping ( $h_s_1, $h_e_1, $h_s_2, $h_e_2, $ext );
	    my ( $tail_o_s, $tail_o_e ) = Overlapping ( $t_s_1, $t_e_1, $t_s_2, $t_e_2, $ext );

	    if ( $head_o_s && $head_o_e && $tail_o_s && $tail_o_e ){
		my $pet_count = $count_1 + $count_2;
		push @clustered_pets, join("-", ( $head_o_s, $head_o_e, $tail_o_s, $tail_o_e, $pet_count ));
		$head_overlap_pets[$i] = 'NA';

		$overlap_flag = 1;

		goto PET;
	    }

	}

	push @clustered_pets, $pet_1;
    }


    while ( $overlap_flag ){

	@head_overlap_pets = @clustered_pets;
	undef @clustered_pets;
	undef $overlap_flag;

	goto PET;

    }

# print clustered pets
    for my $clustered_pet ( @clustered_pets ){
	my @p = split(/-/, $clustered_pet);

	my ( $head_anchor_start, $head_anchor_end ) = ( $p[0], $p[1] );
	my ( $tail_anchor_start, $tail_anchor_end ) = ( $p[2], $p[3] );
	my $cluster_count = $p[4];

## only print clusters
	if ( $cluster_count >= 2 ){
	    print $fh_out join( "\t", ( $chr_head, $head_anchor_start, $head_anchor_end, $chr_tail, $tail_anchor_start, $tail_anchor_end, $cluster_count ) ), "\n";
	}
	else {
	    next;
	}
    }

    return 0;
}



## fix bug in this function
sub Overlapping {
    my ( $s_1, $e_1, $s_2, $e_2, $ext ) = @_;

    my $ss_1 = $s_1 - $ext;
    my $ee_1 = $e_1 + $ext;

    my $ss_2 = $s_2 - $ext;
    my $ee_2 = $e_2 + $ext;

    my @cor = ( $s_1, $e_1, $s_2, $e_2 );
    @cor = sort { $a <=> $b } @cor;

# sorting condition order to improve speed

    if ( ($ss_2 >= $ss_1 && $ss_2 <= $ee_1) && $ee_2 >= $ee_1 ){
	return ( $cor[0], $cor[-1] );
    }
    elsif ( $ss_2 >= $ss_1 && $ee_2 <= $ee_1 ){
	return ( $cor[0], $cor[-1] );
    }
    elsif ( $ss_2 > $ee_1 ){
	return ( 0, 0 );
    }
    ### not happend
    elsif ( $ss_2 <= $ss_1 && $ee_2 >= $ee_1 ){
	return ( $cor[0], $cor[-1] );
    } 
    elsif ( $ee_2 < $ss_1 ){
	return ( 0, 0 );
    }
    elsif ( $ss_2 < $ss_1 && ($ee_2 >= $ss_1 && $ee_2 <= $ee_1) ){
	return ( $cor[0], $cor[-1] );
    } 
    else {
	print STDERR 'Errow in overlapping: ', join ( "\t", ($s_1, $e_1, $s_2, $e_2) ), "\n";
	exit 1;
    }
}



## fix a bug for filtering self-ligation
## Apr 17th, 2017

sub Parsing_Indexed_Bedpe {
    my ( $input, $chr_head, $chr_tail, $self_ligation_cutoff, $large_span_cutoff, $CMD_pairix ) = @_;

    my %pets;
    
    open( my $fh_in, "$CMD_pairix $input \'$chr_head\|$chr_tail\' | " ) or die $!;

    while( my $line = <$fh_in> ){
	chomp $line;
	my @lines = split( /\t/, $line );
	
## exclude read without linker
## only include the reads with linker
	next if $lines[6] =~ /^LinkerNo\:/;

## inter-chr
	if( $chr_head ne $chr_tail ){
	    push @{ $pets{ $lines[1] }{ $lines[2] } }, [ $lines[4], $lines[5] ];
	}
## intra-chr
	else{
	    my $span = $lines[5] - $lines[1];

	    if( $span >= $self_ligation_cutoff && $span <= $large_span_cutoff ){
		push @{ $pets{ $lines[1] }{ $lines[2] } }, [ $lines[4], $lines[5] ];
	    }
	}
    }
    
    close $fh_in;
    return \%pets;
}



sub Create_Chr_Pairs {
    my $chr_names       = shift;
    my $clustering_type = shift;

    my @chr_pairs;

    if( $clustering_type eq 'all' ){
	while( my $chr_1 = shift @{ $chr_names } ){
	    for my $chr_2 ( $chr_1, @{ $chr_names } ){
		push @chr_pairs, [$chr_1, $chr_2];
	    }
	}
    }
    elsif( $clustering_type eq 'intra' ){
	while( my $chr_1 = shift @{ $chr_names } ){
		push @chr_pairs, [$chr_1, $chr_1];
	}
    }
    elsif( $clustering_type eq 'inter' ){
	while( my $chr_1 = shift @{ $chr_names } ){
	    for my $chr_2 ( @{ $chr_names } ){
		push @chr_pairs, [$chr_1, $chr_2];
	    }
	}
    }
    else {
	print STDERR "Unknow parameter $clustering_type\n";
	print STDERR "Accept parameter is all, intra or inter\n";
	exit 1;
    }
    return \@chr_pairs;
}


sub Merge_Files {
    my ( $files, $output ) = @_;

    open( my $fh_out, ">", $output ) or die $!;

    for my $f ( @{ $files } ){
	open( my $fh_in, $f ) or die $!;
	while( my $line = <$fh_in> ){
	    print $fh_out $line;
	}
	close $fh_in;
	unlink( $f );
    }
    close $fh_out;
    return 0;
}


sub Random_String {
    my @chars = ( 'A'..'Z', 'a'..'z' );
    my $string;
    $string .= $chars[rand @chars] for 1..8;

    return $string;
}
