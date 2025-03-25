#!/usr/bin/env perl

chomp($date=`date`);

# inserting the date into ccx_2.15.c

@ARGV="src/ccx_2.15.c";
$^I=".old";
while(<>){
    s/You are using an executable made on.*/You are using an executable made on $date\\n");/g;
    print;
}

@ARGV="src/frd.c";
$^I=".old";
while(<>){
    s/COMPILETIME.*/COMPILETIME       $date                    \\n\",p1);/g;
    print;
}

system "rm -f ccx_2.15.c.old";
system "rm -f frd.c.old";
