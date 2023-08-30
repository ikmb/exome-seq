#!/usr/bin/env ruby
# == NAME
# script_skeleton.rb
#
# == USAGE
# ./this_script.rb [ -h | --help ]
#[ -i | --infile ] |[ -o | --outfile ] |
# == DESCRIPTION
# A skeleton script for Ruby
#
# == OPTIONS
# -h,--help Show help
# -i,--infile=INFILE input file
# -o,--outfile=OUTFILE : output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#  Marc Hoeppner, mphoeppner@gmail.com

require 'optparse'
require 'ostruct'


### Define modules and classes here


### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "A script description here"
opts.separator ""
opts.on("-f","--fasta", "=FASTA","Fasta file") {|argument| options.fasta = argument }
opts.on("-b","--base", "=BASE","Base name") {|argument| options.base = argument }
opts.on("-r","--report", "=REPORT","Assembly report") {|argument| options.report = argument }
opts.on("-o","--outfile", "=OUTFILE","Output file") {|argument| options.outfile = argument }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse!

file = File.open(options.fasta,"r")

base = options.base

out = File.new("#{base}.fasta","w")
l = File.new("#{base}_lookup.txt","w")

# are we building a non_alt reference?
base.include?("_alt") ? alt = true : alt = false

report = IO.readlines(options.report)

lookup = {}

report.each do |line|

    next if line.match(/^#.*/)

    e = line.strip.split("\t").collect{|m| m.strip }

    acc = e[6]
    ucsc = e[-1]

    if ucsc == "na"
        t = e[1]
        a = e[4]
        c = e[2]
        base = "chr#{c}_#{a.upcase.gsub('.','v')}"
        t == "novel-patch" ? ucsc = "#{base}_alt" : ucsc = "#{base}_fix"
        warn "Renamed na to #{ucsc} (#{t})"
    end

    l.puts "#{acc}\t#{ucsc}"

    lookup[acc] = ucsc

end

l.close

skip = false

while (line = file.gets)

    line.strip!

    if line.match(/^>.*/)
        this_acc = line.strip.split(" ")[0].gsub(">", "")
        if lookup.has_key?(this_acc)
            acc = lookup[this_acc]

            # if we are building a non-alt ref
            if alt
                # skip all entries that are either ALT or FIX patches
                if acc.include?("_alt") || acc.include?("_fix") 
                    skip = true 
                else
                    skip = false
                end
            end

            out.puts ">#{acc}" unless skip
        
        else
                abort line
        end

    else

        out.puts line unless skip

    end

end

file.close
out.close
