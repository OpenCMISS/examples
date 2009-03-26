#!/usr/bin/env perl

use warnings;
use strict;
use Digest::MD5;

sub md5sum{
  my $file = shift;
  my $digest = "";
  eval{
    open(FILE, $file) or die "Can't find file $file\n";
    my $ctx = Digest::MD5->new;
    $ctx->addfile(*FILE);
    $digest = $ctx->hexdigest;
    close(FILE);
  };
  if($@){
    print $@;
    return "";
  }
  return $digest;
}

sub usage{
  print "usage: ./testscript.pl path-to_executable\n";
  exit 1;
}

if($#ARGV + 1 != 1){
  usage();
 }

my $cmd = $ARGV[0];
system $cmd;

my $actualfile = "LaplaceExample";
my $expectedfile = "ExpectedLaplaceExample";
my $md5actual =  md5sum($actualfile);
my $md5expected =  md5sum($expectedfile);
if($md5actual eq $md5expected){
  print "Analytic Laplace Example Testcase3 - Output file check is successfully completed.\n";
}else{
  print "ERROR:The output file is not as expected!\n";
  exit 1;
}
exit 0;
