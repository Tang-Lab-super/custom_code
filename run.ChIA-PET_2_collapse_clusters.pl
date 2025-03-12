#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw ( max sum );
use threads;
use threads::shared;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use GenomeSize;

## updated Jun 3rd 2019

## updated Apr 7th, 2017
## change the anchor extenion to 0 bp

## updated Apr 24th, 2017
## include loops with broad anchors
## ERROR Anchor orientation:
## chr14   107026461       107246095       chr14   107096761       107271603       9374

## fixed the chromosome sorting problem
## sub Sort_by_Anchor_Position is abandoned

my $input_clusters        = shift;
my $output_mergedClusters = shift;
my $ref_genome            = shift;
my $cpu_num               = shift;

####################### default extension for clustering ###################################
my $anchor_ext = 0;

###################### reference genome chromosome index ##################################
my $genome    = GenomeSize->new( genome => $ref_genome );
my $chr_index = $genome->Get_chr_index();
my @chr_names = sort { $chr_index->{ $a } <=> $chr_index->{ $b } } keys %{ $chr_index };


################################ generate random prefix ###################################
my $prefix = Random_String();

###########################################################################################


################################### mergeing clusters #####################################
my $anchors_pos    = Parsing_PetCluster ( $input_clusters );
my $cluster        = Parsing_PetCluster_File ( $input_clusters );
my $anchor_cluster = Anchor_Position_Extend ( $anchors_pos, $anchor_ext, $chr_index );

print STDERR "Merging clusters ...\n";
my $output_mergedCluster_files = Run_Merging_Clusters_Multi_Threads ( \@chr_names, $anchors_pos, $cluster, $anchor_cluster, $cpu_num, $prefix );

print STDERR "Merging output files ...\n";
&Merge_Files ( $output_mergedCluster_files, $output_mergedClusters );





####################################### subrotinates #########################################

sub Run_Merging_Clusters_Multi_Threads {
    my ( $chr_names_ref, $anchors_pos, $cluster, $anchor_cluster, $cpu_num, $prefix ) = @_;

    my @output_clusters :shared;
    my $output_clusters_ref = \@output_clusters;

    my @chr_names = @{ $chr_names_ref };

## require all numbers of CPU
    if( scalar @chr_names >= $cpu_num ){
	while( (scalar @chr_names - $cpu_num) >= 0 ){
	    my @threads;
	    my @sub_group = splice( @chr_names, 0, $cpu_num );

	    foreach my $chr ( @sub_group ){
		my $thr = threads->create( 'Run_Merging_Cluster', 
					   ( $chr, $anchors_pos, $cluster, $anchor_cluster, $output_clusters_ref, $prefix ) );
		push @threads, $thr;
	    }
	    foreach ( @threads ){
		my $tid = $_->tid();
		print STDERR "Thread collected: $tid\n";

		$_->join();
		
		if ( @chr_names ){
		    my $chr_add = shift @chr_names;
		    my $thr_add = threads->create( 'Run_Merging_Cluster',
						   ( $chr_add, $anchors_pos, $cluster, $anchor_cluster, $output_clusters_ref, $prefix ) );
		    push @threads, $thr_add;
		}
	    }
	}
    }

## needs few numbers of CPU
    else {
	my @threads;
	for my $chr ( @chr_names ){
	    my $thr = threads->create( 'Run_Merging_Cluster', 
				       ( $chr, $anchors_pos, $cluster, $anchor_cluster, $output_clusters_ref, $prefix ) );
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


sub Run_Merging_Cluster {
    my ( $chr, $anchors_pos, $cluster, $anchor_cluster, $output_clusters_ref, $prefix ) = @_;
    
    
    my @mergedAnchor_regions = keys %{ $anchor_cluster->{ $chr } };
#    print Dumper( @mergedAnchor_regions );

    @mergedAnchor_regions = map { $_->[0] }
                           sort { $a->[-2] <=> $b->[-2] || $a->[-1] <=> $b->[-1] }
                           map { [$_, split(/-/, $_)] } @mergedAnchor_regions;

    my $output = $prefix.'.'.$chr.'.merged.clusters.txt';
    push @{ $output_clusters_ref }, $output;

    open( my $fh_out, ">", $output ) or die $!;

    while ( my $target_mAnchor = shift @mergedAnchor_regions ){

	for my $query_mAnchor ( @mergedAnchor_regions ) {
	    my ( $overlapped_clusters, $total_PetCount ) = Overlapping_PetCluster_v2 ( $chr, $target_mAnchor, $query_mAnchor, $anchor_cluster, $cluster );

	    if ( $total_PetCount ){
		my $mAnchor_t_pos = $chr.':'.$target_mAnchor;
		my $mAnchor_q_pos = $chr.':'.$query_mAnchor;

		my $overlapped_PetCluster_infor = join( "|", @{ $overlapped_clusters } );

		print $fh_out join( "\t", ( $mAnchor_t_pos, $mAnchor_q_pos, $total_PetCount, $overlapped_PetCluster_infor )), "\n";

	    }
	}
	
    }
    close $fh_out;

    return $output_clusters_ref;
}


sub Overlapping_PetCluster_v2 {
    my ( $chr, $target_mAnchor, $query_mAnchor, $anchor_cluster, $cluster ) = @_;

    my @target_anchors = @{ $anchor_cluster->{ $chr }{ $target_mAnchor } };
    my @query_anchors  = @{ $anchor_cluster->{ $chr }{ $query_mAnchor }  };

    my @overlapped_PetClusters;
    my $total_PetCounts;

    for my $t ( @target_anchors ) {
	for my $q ( @query_anchors ) {
	    my $petCount = $cluster->{ $chr }{ $t }{ $q };

	    if( $petCount ) {
		push @overlapped_PetClusters, join( ",", ( $chr.':'.$t, $chr.':'.$q, $petCount ) );
		$total_PetCounts += $petCount;
	    }
	}
    }

    return ( \@overlapped_PetClusters, $total_PetCounts );
}




sub Anchor_Position_Extend {
    my ( $anchors_pos, $span, $chr_index ) = @_;

    my ( $chr, $start, $end, %extend_anchor, @connected_anchors, @PetCount );

    for my $anchor ( map { $_->[-1] }
		     sort { $chr_index->{$a->[0]} <=> $chr_index->{$b->[0]} || $a->[1] <=> $b->[1] }
		     map { [split(/[:-]/, $_), $_] }
		     keys %{ $anchors_pos } ) {

	my @anchor_infor = split(/[:-]/, $anchor);

	my $anchor_chr = $anchor_infor[0];
	my $anchor_s   = $anchor_infor[1];
	my $anchor_e   = $anchor_infor[2];
	my $anchor_PetCount = $anchors_pos->{$anchor};

	unless ( $chr && $start && $end && @connected_anchors ){
	    $chr   = $anchor_chr;
	    $start = $anchor_s;
	    $end   = $anchor_e;
	    push @connected_anchors, $anchor;
	    push @PetCount, $anchor_PetCount;
	    next;
	}


	if ( $chr eq $anchor_chr && ( $end + $span ) >= $anchor_s && $anchor_e >= $end ){
	    push @connected_anchors, $anchor;
	    push @PetCount, $anchor_PetCount;
	    $end = $anchor_e;

	} elsif ( $chr eq $anchor_chr && ( $end + $span ) >= $anchor_s && $anchor_e <= $end ){
	    push @connected_anchors, $anchor;
	    push @PetCount, $anchor_PetCount;
	
	} elsif ( $chr eq $anchor_chr && ( $end + $span ) < $anchor_s ){
	    $extend_anchor{$chr . ':' . $start . '-' . $end}{'anchors'}  = [ @connected_anchors ];
	    $extend_anchor{$chr . ':' . $start . '-' . $end}{'PetCount'} = [ @PetCount ];

	    $start = $anchor_s;
	    $end   = $anchor_e;
	    undef @connected_anchors;
	    undef @PetCount;
	    push @connected_anchors, $anchor;
	    push @PetCount, $anchor_PetCount;

	} elsif ( $chr ne $anchor_chr && $start && $end && @connected_anchors ){
	    $extend_anchor{$chr . ':' . $start . '-' . $end}{'anchors'}  = [ @connected_anchors ];
	    $extend_anchor{$chr . ':' . $start . '-' . $end}{'PetCount'} = [ @PetCount ];

	    $chr   = $anchor_chr;
	    $start = $anchor_s;
	    $end   = $anchor_e;
	    undef @connected_anchors;
	    undef @PetCount;
	    push @connected_anchors, $anchor;
	    push @PetCount, $anchor_PetCount;
	} else {
 	    print STDERR 'Unknown case: ', join("\t", ( $anchor_chr, $anchor_s, $anchor_e )), "\n";
	}
    }

    if ( @connected_anchors && @PetCount ){
	    $extend_anchor{$chr . ':' . $start . '-' . $end}{'anchors'}  = [ @connected_anchors ];
	    $extend_anchor{$chr . ':' . $start . '-' . $end}{'PetCount'} = [ @PetCount ];
	    undef @connected_anchors;
	    undef @PetCount;
	    undef $chr;
	    undef $start;
	    undef $end;
    }

    my %extend_anchor_infor;
    
    for my $merged_anchor ( 			    
	map { $_->[-1] }
	sort { $chr_index->{$a->[0]} <=> $chr_index->{$b->[0]} || $a->[1] <=> $b->[1] }
	map { [split(/[:-]/, $_), $_] }
	keys %extend_anchor ){

	my ( $chr, $anchor_pos ) = split(/\:/, $merged_anchor);

	my @cluster_anchors = @{ $extend_anchor{ $merged_anchor }{'anchors' } };
	@cluster_anchors = map{ my @a = split(/\:/, $_); $a[1] } @cluster_anchors;

	$extend_anchor_infor{ $chr }{ $anchor_pos } = \@cluster_anchors;
    }

    return \%extend_anchor_infor;
}



sub Parsing_PetCluster_File {
    my $file = shift;

    my %cluster;

    open ( my $fh_in, $file ) or die $!;
    while ( my $line = <$fh_in> ){
	chomp $line;
	next if $line =~ /^chrM/;
	my @lines = split( /\t/, $line );

## exclude those loops with broad anchors
	if ( $lines[2] > $lines[4] ){
	    print STDERR "Broad anchors (abundant interaction region, but skip this region for merging clusters): \n";
	    print STDERR $line, "\n";
	    next;

	}

	$cluster{ $lines[0] }{ $lines[1].'-'.$lines[2] }{ $lines[4].'-'.$lines[5] } = $lines[6];
    }
    close $fh_in;
    return \%cluster;
}



sub Parsing_PetCluster {
    my $in = shift;
    my %anchors;

    open(my $fh_in, $in) or die "Can't open input cluster file: $!";

    while (my $line = <$fh_in>){
	chomp $line;
	my @lines = split(/\t/, $line);

	next if ( $lines[0] =~ /[C|c]hr[M|m]/ || $lines[3] =~ /[C|c]hr[M|m]/ );

## head anchor position (chr1:3657531-3658374) and PET count
	$anchors{ $lines[0].':'.$lines[1].'-'.$lines[2] } = $lines[6];
## tail anchor position and PET count
	$anchors{ $lines[3].':'.$lines[4].'-'.$lines[5] } = $lines[6];
    }
    close $fh_in;
    return \%anchors;
}



sub Sort_by_Anchor_Position {
    my @a = split(/[:-]/, $a);
    my @b = split(/[:-]/, $b);

    $a[0] =~ /chr(.+)/;
    my $a_chr_index = $1;

    $b[0] =~ /chr(.+)/;
    my $b_chr_index = $1;

    if ( $a_chr_index eq 'X' or $a_chr_index eq 'x' ){
	push @a, 98;
    } elsif ( $a_chr_index eq 'Y' or $a_chr_index eq 'y' ){
	push @a, 99;
    } elsif ( $a_chr_index eq 'M' or $a_chr_index eq 'm' ){
	push @a, 100;
    } elsif ( $a_chr_index > 0 ){
	push @a, $a_chr_index;
    } else {
	print STDERR "Something wrong in chr index: $a_chr_index\n";
	exit 1;
    }

    if ( $b_chr_index eq 'X' or $b_chr_index eq 'x' ){
	push @b, 98;
    } elsif ( $b_chr_index eq 'Y' or $b_chr_index eq 'y' ){
	push @b, 99;
    } elsif ( $b_chr_index eq 'M' or $b_chr_index eq 'm' ){
	push @b, 100;
    } elsif ( $b_chr_index > 0 ){
	push @b, $b_chr_index;
    } else {
	print STDERR "Something wrong in chr index: $b_chr_index\n";
	exit 1;
    }

    $a[-1] <=> $b[-1]
	   ||
    $a[1] <=> $b[1]
           ||
    $a[2] <=> $b[2]
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
