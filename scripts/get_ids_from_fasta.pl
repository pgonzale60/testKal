#!/usr/bin/env perl

use warnings;
use strict;

## Take a list of IDs, obtain fasta sequences from bigger file

use strict;
my $rem_polyA = 1;   # remove run of As at end of sequences? Useful for RefSeq 3'UTRs
my $min_len   = 50;
my $sort_mode = 0;

my $usage = "Usage: $0 id_list fasta_file [min_len] [rem_polyA:1|0] [sort_by_length:0|1]";

unless (scalar @ARGV >= 2) {
   die "$usage\n";
}
my ($ids,$fasta,$len,$polyA,$sort_by_length) = @ARGV;
unless (-s $ids and -s $fasta) {
   die "$usage\n";
}
$min_len = $len if $len;
$rem_polyA = $polyA if $polyA;
$sort_mode = $sort_by_length if $sort_by_length;


## Open ids
open (IDS, "gzip -dcf $ids |") || die "Can't gzip -dcf $ids!\n";
my %ids;
my @order;
if ($sort_mode) {
   while (<IDS>) {
      next unless /^(\S+)/;
      $ids{$1} = 1;
   }
} else {
   while (<IDS>) {
      next unless /^(\S+)/;
      $ids{$1} = 1;
      push @order, $1;
   }
}
close IDS;

## Read fasta sequences, print out only those from ids
open (FASTA, "gzip -dcf $fasta |") || die "Can't gzip -dcf $fasta!\n";
$/ = ">";
<FASTA>;
my %all_seqs;
while (<FASTA>) {
   chomp;
   next unless s/^(\S+)\.?.*\n//;
   my $id = $1;
   s/\s+//g;
   my $seq = uc $_;
   if ($rem_polyA) {
      $seq =~ s/A{2,}$//;
   }
   next unless length($seq) >= $min_len;
   foreach my $id (split /[\|,]/, $id) {
      next unless $ids{$id};
      $all_seqs{$id} = $seq;
   }
}
$/ = "\n";
close FASTA;

if ($sort_mode){
   foreach my $id (sort {length($all_seqs{$a}) <=> length($all_seqs{$b})} keys %all_seqs) {
      print ">$id\n$all_seqs{$id}\n";
   }
} else {
   my @unique = uniq(@order);
   foreach my $id (@unique) {
      if ($all_seqs{$id}) {
         print ">$id\n$all_seqs{$id}\n";
      }
   }
}

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
