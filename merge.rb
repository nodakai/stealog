#! /usr/bin/env ruby

h = Hash.new() {|h,k| h[k] = [] }
id = nil
ATMARK = '@'[0]
ARGF.readlines.each {|l|
    if l[0] == ATMARK
        l =~ %r[tabelog/(.*)_dtlrvwlst[^,]*,(.*)\Z]
        id = $1 + ',' + $2
    elsif l =~ /^\.\./
        raise "Unexpected line \"#{l}\""
    else
        h[id] << l
    end
}

STDERR.puts h.size
# p h.values.max {|a,b| a.size <=> b.size }

h.each {|id,revArr|
    if revArr.size > 1
        puts "@" + id
        revArr.each {|l| puts l }
    end
}
