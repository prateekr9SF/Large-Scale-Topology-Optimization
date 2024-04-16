#!/usr/bin/env perl

chomp($date=`date`);

# inserting the date into main.c

@ARGV="main.c";
$^I=".old";
while(<>){
    s/You are using an executable made on.*/You are using an executable made on $date\\n");/g;
    print;
}

@ARGV="frd.c";
$^I=".old";
while(<>){
    s/COMPILETIME.*/COMPILETIME       $date                    \\n\",p1);/g;
    print;
}

system "rm -f main.c.old";
system "rm -f frd.c.old";
