#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin qw($Bin);
#use lib "$Bin/../lib";
use lib "/data/public/software/CHIAPET_PIPELINE/lib/";
use GenomeSize;
use BSearch qw( Bsearch_num_range );

my $clusters_in    = shift; # in txt formate
my $peaks_in       = shift; # in bed formate, summit is the middle of peak region
my $ref_genome     = shift;
my $output_cluster = shift;
my $output_peak    = shift;

my $ext = 1000;

my $genome    = GenomeSize->new( genome => $ref_genome );
my $chr_index = $genome->Get_chr_index();


#my $clusters  = Parsing_Clusters ( $clusters_in );

my $clusters  = Parsing_Merged_Clusters ( $clusters_in );
my $peaks     = Parsing_Peaks           ( $peaks_in    );

my %overlapped_peaks;

open ( my $fh_out_clusters, ">", $output_cluster ) or die $!;
open ( my $fh_out_peaks,    ">", $output_peak    ) or die $!;

print $fh_out_clusters join( "\t", ( 'HeadAnchor_pos', 'TailAnchor_pos', 'PetCount', 'Type_support', 'HeadAnchor_peak_summit', 'TailAnchor_peak_summit' ) ), "\n";

for my $chr ( sort { $chr_index->{ $a } <=> $chr_index->{ $b } } keys %{ $chr_index } ){

    my $peak_boundary = Get_Peak_Boundary ( $peaks, $chr );

    for my $headAnchor ( map { $_->[0] }
			 sort { $a->[-2] <=> $b->[-2] || $a->[-1] <=> $b->[-1] }
			 map { [ $_, split(/-/, $_) ] }
			 keys %{ $clusters->{ $chr } } ){

	for my $tailAnchor ( map { $_->[0] }
			     sort { $a->[-2] <=> $b->[-2] || $a->[-1] <=> $b->[-1] }
			     map { [ $_, split(/-/, $_) ] }
			     keys %{ $clusters->{ $chr }{ $headAnchor } } ){

	    my $freq = $clusters->{ $chr }{ $headAnchor }{ $tailAnchor };

	    my ( $headAnchor_s, $headAnchor_e ) = split( /-/, $headAnchor );
	    my ( $tailAnchor_s, $tailAnchor_e ) = split( /-/, $tailAnchor );

	    my $headAnchor_s_inx = Bsearch_num_range ( ( $headAnchor_s - $ext ), $peak_boundary, 's' );
	    my $headAnchor_e_inx = Bsearch_num_range ( ( $headAnchor_e + $ext ), $peak_boundary, 'e' );

	    my $tailAnchor_s_inx = Bsearch_num_range ( ( $tailAnchor_s - $ext ), $peak_boundary, 's' );
	    my $tailAnchor_e_inx = Bsearch_num_range ( ( $tailAnchor_e + $ext ), $peak_boundary, 'e' );

	    my $headAnchor_peaks = Get_Overlapped_Peaks ( $headAnchor_s_inx, $headAnchor_e_inx, $peak_boundary, $peaks, $chr );
	    my $tailAnchor_peaks = Get_Overlapped_Peaks ( $tailAnchor_s_inx, $tailAnchor_e_inx, $peak_boundary, $peaks, $chr );

	    my $headAnchor_peaks_infor;
	    my $tailAnchor_peaks_infor;

	    if ( $headAnchor_peaks eq 'NA' ){
		$headAnchor_peaks_infor = 'NA';
	    } else {
		$headAnchor_peaks_infor = join( "\|", @{ $headAnchor_peaks } );
	    }

	    if ( $tailAnchor_peaks eq 'NA' ){
		$tailAnchor_peaks_infor = 'NA';
	    } else {
		$tailAnchor_peaks_infor = join( "\|", @{ $tailAnchor_peaks } );
	    }

	    for my $anchor_peaks ( $headAnchor_peaks, $tailAnchor_peaks ){
		next if $anchor_peaks eq 'NA';

		for my $p ( @{ $anchor_peaks } ){
		    my ( $p_chr, $p_s, $p_e, $p_freq ) = split( /[\:\-\,]/, $p );
		    $overlapped_peaks{ $p_chr }{ $p_s.'-'.$p_e } = $p_freq;
		}
	    }


	    my $loop_type = Loop_with_Peak_Support( $headAnchor_peaks_infor, $tailAnchor_peaks_infor );


	    print $fh_out_clusters join( "\t", ( $chr.':'.$headAnchor, $chr.':'.$tailAnchor, $freq, $loop_type, $headAnchor_peaks_infor, $tailAnchor_peaks_infor ) ), "\n";

	} #tailAnchor

    } #headAnchor

}


for my $chr ( sort { $chr_index->{ $a } <=> $chr_index->{ $b } } keys %{ $chr_index } ){
    for my $pos ( map { $_->[0] }
		  sort { $a->[-2] <=> $b->[-2] || $a->[-1] <=> $b->[-1] }
		  map { [$_, split(/-/, $_)] }
		  keys %{ $peaks->{ $chr } } ){

	my $peak_freq = $peaks->{ $chr }{ $pos };

	if ( $overlapped_peaks{ $chr }{ $pos } ){
	    print $fh_out_peaks join( "\t", ( $chr.':'.$pos, $peak_freq, 'Overlapped' ) ), "\n";
	} else {
	    print $fh_out_peaks join( "\t", ( $chr.':'.$pos, $peak_freq, 'None' ) ), "\n";
	}
    } 
}


close $fh_out_clusters;
close $fh_out_peaks;



sub Get_Overlapped_Peaks {
    my ( $s_inx, $e_inx, $list, $peaks, $chr ) = @_;

    if ( $s_inx > $e_inx ){
	print STDERR 'Wrong index boundary to fectch boundary from list:', 
	join( "\t", ( $s_inx, $e_inx ) ), "\n";
	exit 1;
    }

    elsif ( $s_inx == $e_inx ) {
	return 'NA';
    } 
    else {
	my ( $overlapped_s_inx, $overlapped_e_inx );

	if ( $s_inx % 2 == 0 ){
	    $overlapped_s_inx = $s_inx;
	} else {
	    $overlapped_s_inx = $s_inx + 1;
	}

	if ( $e_inx % 2 == 0 ){
	    $overlapped_e_inx = $e_inx - 1;
	} else {
	    $overlapped_e_inx = $e_inx;
	}

	my @overlapped_peaks;
	for ( my $i = $overlapped_s_inx; $i <= $overlapped_e_inx; $i += 2 ){
	    my $peak_s = $list->[$i];
	    my $peak_e = $list->[$i+1];

	    my $peak_freq = $peaks->{ $chr }{ $peak_s.'-'.$peak_e };
	    push @overlapped_peaks, $chr.':'.$peak_s.'-'.$peak_e.','.$peak_freq;

	}
	if ( @overlapped_peaks ){
	    return \@overlapped_peaks;
	} else {
	    return 'NA';
	}
    }
}



sub Binary_Search_Index {
    my $list = shift;
    my $site = shift;
    my $flag = shift;

    my $last_index = scalar @{ $list } - 1;

    my ( $low, $high ) = ( 0, $last_index );

    if ( $site < $list->[0] ){
	return 0;
    } elsif ( $site > $list->[-1] ){
	return $last_index;
    }

    while ( $low < $high ){ #while ( $low <= $high ){

	my $try= int( $low + ( $high - $low ) / 2 );
	
	if ( $list->[$try] == $site ){
	    if ( $flag eq 's' ){
		return $try - 1;
	    } elsif ( $flag eq 'e' ){
		return $try + 1;
	    }
#            return $try;
        }

	elsif ($list->[$try+1] == $site){
	    if ( $flag eq 's' ){
		return $try + 1 - 1;
	    } elsif ( $flag eq 'e' ){
		return $try + 1 + 1;
	    }
#	    return $try + 1;
	}

	elsif ($list->[$try-1] == $site){
	    if ( $flag eq 's' ){
		return $try - 1 - 1;
	    } elsif ( $flag eq 'e' ){
		return $try - 1 + 1;
	    }
#	    return $try - 1;
	}

	elsif ($list->[$try] > $site && $list->[$try-1] < $site){
	    if ( $flag eq 's' ){
		return $try - 1;
	    } elsif ( $flag eq 'e' ){
		return $try - 1 + 1;
	    }
#	    return $try - 1;
	}

	elsif( $list->[$try] < $site && $list->[$try+1] > $site){
	    if ( $flag eq 's' ){
		return $try;
	    } elsif ( $flag eq 'e' ){
		return $try + 1;
	    }
#	    return $try;
	} 

#	elsif ( $list->[$try] < $site && $list->[$try + 1] < $site ) {
#            $low = $try + 1;
#            next;
#        }
	elsif ( $list->[$try] < $site  ) {
            $low = $try;
            next;
        }

#	elsif( $list->[$try] > $site && $list->[$try - 1] > $site){
#            $high = $try -1;
#            next;
#        }
	elsif( $list->[$try] > $site ){
            $high = $try;
            next;
        }
    }
}




sub Get_Peak_Boundary {
    my $peaks = shift;
    my $chr   = shift;

    print STDERR 'Calculate sorted list of peak boundaries ...', "\n";
    my @boundary_sorted;

    for my $peak_pos ( map { $_->[0] }
		       sort { $a->[-2] <=> $b->[-2] || $a->[-1] <=> $b->[-1] }
		       map { [$_, split(/-/, $_)] }
		       keys %{ $peaks->{ $chr } } ){

	my ( $peak_s, $peak_e ) = split( /-/, $peak_pos );
	push @boundary_sorted, ( $peak_s, $peak_e );

    }


    if ( $boundary_sorted[0] ){
	return \@boundary_sorted;
    } else {
	push @boundary_sorted, 0;
	return \@boundary_sorted;
    }

}


sub Parsing_Peaks {
    my $in = shift;

    print STDERR 'Parsing peaks ...', "\n";

    open ( my $fh_in, $in ) or die $!;

    my %peaks;

    while ( my $line = <$fh_in> ){
	chomp $line;
	next if $line =~ /^\#/;
	next if $line =~ /^$/;
	next if $line =~ /^chr[\t|\s]+start/;

	my @lines = split( /[\t|\s]/, $line );

	next if ( $lines[0] eq 'chrM' );

	my $peak_chr = $lines[0];
	my $peak_pos = $lines[1].'-'.$lines[2];
	my $summit   = $lines[1] + int(($lines[2]-$lines[1])/2);

	$peaks{ $peak_chr }{ $peak_pos } = $summit;
    }
    close $fh_in;
    return \%peaks;
}




sub Loop_with_Peak_Support {
    my $head_peak = shift;
    my $tail_peak = shift;

    my $type;

    if( $head_peak =~ /^chr/ && $tail_peak =~ /^chr/ ){
	$type = 'Both';
    }
    elsif ( ($head_peak =~ /^chr/ && $tail_peak eq 'NA') || ($head_peak eq 'NA' && $tail_peak =~ /^chr/) ){
	$type = 'One';
    }
    elsif ( $head_peak eq 'NA' && $tail_peak eq 'NA' ){
	$type = 'None';
    }
    else {
	print STDERR "Unknow type of peak support!\n";
	print STDERR "\t$head_peak\t$tail_peak\n";
	exit 1;
    }

    return $type;
}


sub Parsing_Merged_Clusters {
    my $input = shift;

    my $fh_in = Create_Input_FH ( $input );

    my %clusters;

    while ( my $line = <$fh_in> ){
	chomp $line;
	my @lines = split( /\t/, $line );

	my $freq = $lines[2];

	my ( $chr_head, $start_head, $end_head ) = split( /[:-]/, $lines[0] );
	my ( $chr_tail, $start_tail, $end_tail ) = split( /[:-]/, $lines[1] );

## only consider intra-chr interactions
	if( $chr_head eq $chr_tail ){
	    my $headAnchor = $start_head.'-'.$end_head;
	    my $tailAnchor = $start_tail.'-'.$end_tail;
	    $clusters{ $chr_head }{ $headAnchor }{ $tailAnchor } = $freq;	    
	}
	else {
	    next;
	}
    }

    close $fh_in;
    return \%clusters;
}



## create input file handle
sub Create_Input_FH {
    my $input = shift;
    my $fh;

    if( $input =~ /\.gz$/ ){
	open( $fh, "gunzip -c $input |" ) or die $!;
	return $fh;
    }
    else {
	open( $fh, $input ) or die $!;
	return $fh;
    }
}
