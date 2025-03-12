#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use FindBin qw($Bin);
#use lib "$Bin/../lib";
use lib "/data/public/software/CHIAPET_PIPELINE/lib/";
use BSearch qw( Bsearch_num_range );
use Overlapping qw(OverlappingTest);

## change output of the script
## May 30, 2017

my $input_clusters  = shift; 
## merged cluster in this format: chr1:1875012-1876376    chr1:1975716-1980086    8       chr1:1875012-1876190,chr1:1977865-1980086,5|chr1:1875515-1876376,chr1:1975716-1976556,3
my $input_bedgraph  = shift; ## bedgraph for coverage
my $output_clusters = shift;


my $COVERAGE_BASELINE = 1;

print STDERR "Reading bedgraph ...\n";
my $coverage = Parsing_Bedgraph ( $input_bedgraph, $COVERAGE_BASELINE );

print STDERR "Reading clusters ...\n";
my ( $cluster, $anchor ) = Parsing_Cluster ( $input_clusters );

print STDERR "Adjusting anchor position ...\n";
my $anchor_adjusted = Adjust_Anchor_Pos ( $anchor, $coverage );


print STDERR "Print cluster with adjust anchors ...\n";

open( my $fh_out, ">", $output_clusters ) or die $!;
print $fh_out join( "\t", ( 'headAnchor_adjust', 'tailAnchor_adjust', 'petCount', 'headAnchor', 'tailAnchor', 'petCount', 'peakSupport', 'headAnchor_peak', 'tailAnchor_peak' ) ), "\n";

for my $head ( keys %{ $cluster } ){
    for my $tail ( keys %{ $cluster->{ $head } } ){
	my $count        = $cluster->{ $head }{ $tail }->[0];
	my $peak_support = $cluster->{ $head }{ $tail }->[1];
	my $peak_head    = $cluster->{ $head }{ $tail }->[2];
	my $peak_tail    = $cluster->{ $head }{ $tail }->[3];
		
	my $head_adjust = $anchor_adjusted->{ $head };
	my $tail_adjust = $anchor_adjusted->{ $tail };

	print $fh_out join( "\t", ( $head_adjust, $tail_adjust, $count, $head, $tail, $count, $peak_support, $peak_head, $peak_tail ) ), "\n";
    }
}
close $fh_out;


sub Adjust_Anchor_Pos {
    my ( $anchor, $coverage ) = @_;

    my %adjust_anchor;

    for my $pos ( keys %{ $anchor } ){

	my $adjust_anchor_coverage = Adjust_Anchor_Pos_Coverage( $pos, $coverage );

	if( $adjust_anchor_coverage ){
	    $adjust_anchor{ $pos } = $adjust_anchor_coverage;
	}
	else {
	    $adjust_anchor{ $pos } = $pos;
	}
    }

    return \%adjust_anchor;
}

sub Adjust_Anchor_Pos_Coverage {
    my ( $pos, $coverage ) = @_;
    
    my ( $chr, $start, $end ) = split( /[:-]/, $pos );
    
    my $list_cov_start = $coverage->{ $chr }{ 'start' };
    my $list_cov_end   = $coverage->{ $chr }{ 'end'   };

    my $start_index_cov = Bsearch_num_range ( $start, $list_cov_start, 's' );
    my $end_index_cov   = Bsearch_num_range ( $end,   $list_cov_end, 'e' );

    my $ext_index = 5;
    my $start_index_cov_ext = $start_index_cov - $ext_index;
    my $end_index_cov_ext   = $end_index_cov   + $ext_index;

    my @overlapped_reg;

## overlapping search
    for my $index ( $start_index_cov_ext..$end_index_cov_ext ){
	my $reg_s = $coverage->{ $chr }{ 'start' }->[$index];
	my $reg_e = $coverage->{ $chr }{ 'end' }->[$index];
	my $reg_v = $coverage->{ $chr }{ 'cov' }->[$index];

	next unless ( $reg_s && $reg_e && $reg_v );
	
	my $overlapping = OverlappingTest( $start, $end, $reg_s, $reg_e );

	if( $overlapping ){
	    push @overlapped_reg, [$reg_s, $reg_e, $reg_v];
	}
    }


    if( $overlapped_reg[0] ){
	my @overlapped_reg_sorted = sort { $b->[-1] <=> $a->[-1] } @overlapped_reg;
	my $adjust_pos = $chr.':'.$overlapped_reg_sorted[0]->[0].'-'.$overlapped_reg_sorted[0]->[1];
	return $adjust_pos;
    }
    else {
	return 0;
    }

}


sub Parsing_Cluster {
    my $input = shift;
    
    my $fh_in = Create_Input_FH( $input );

    my ( %cluster, %anchor );
    while( my $line = <$fh_in> ){
	next unless $line =~ /^chr/;

	chomp $line;
	my @lines = split( /\t/, $line );
	
	my $head_pos     = $lines[0];
	my $tail_pos     = $lines[1];
	my $count        = $lines[2];
	my $peak_support = $lines[3];
	my $peak_head    = $lines[4];
	my $peak_tail    = $lines[5];

	$cluster{ $head_pos }{ $tail_pos } = [$count, $peak_support, $peak_head, $peak_tail];

	$anchor{ $head_pos }++;
	$anchor{ $tail_pos }++;
    }
    close $fh_in;

    return ( \%cluster, \%anchor );
}


sub Parsing_Bedgraph {
    my $input    = shift;
    my $baseline = shift;
    
    my $fh_in = Create_Input_FH( $input );

    my %cov;
    while( my $line = <$fh_in> ){
	chomp $line;
	my @lines = split( /\t/, $line );
	if ( $lines[3] > $baseline ){
	    push @{ $cov{ $lines[0] } }, [$lines[1], $lines[2], $lines[3]];
	}
    }
    close $fh_in;


    my %cov_sorted;
    for my $chr ( keys %cov ){
	my @pos_sorted = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{ $cov{ $chr } };
	push @{ $cov_sorted{ $chr }{'start'} }, map{ $_->[0] } @pos_sorted;
	push @{ $cov_sorted{ $chr }{'end'} }, map{ $_->[1] } @pos_sorted;
	push @{ $cov_sorted{ $chr }{'cov'} }, map{ $_->[2] } @pos_sorted;

	delete $cov{ $chr };
    }

    return \%cov_sorted;
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
