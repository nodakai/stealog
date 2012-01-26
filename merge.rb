#! /usr/bin/env ruby

h = Hash.new() {|h,k| h[k] = [] }
id = nil
ARGF.readlines.each {|l|
    case l
    when %r[tabelog/(.*)_dtlrvwlst]
        id = $1
    when /^\.\./
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
