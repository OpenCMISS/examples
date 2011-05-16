# extracts component 4-6 from nodes (deformed coords) to fix broken .exnode files from fluid io
Dir['*.exnode'].each{|f|
  File.open(f.sub('S_TIMESTP','OUT'),'w'){|fp|
    fp << (File.read(f)+'END').sub(/Fields=\s*\d+/,'Fields= 1').sub(/\s*2\).*?Node/m,"\nNode").gsub(/(Node:\s*\d+\s*)(.*?)(?=Node|END)/m){|m|  $1 + $2.split(/\s+/).map{|l| l.to_f}[3..5].join("\n       ")+"\n"}[0..-4]
  } unless f['OUT']
}
Dir['*.exelem'].each{|f| 
  File.open(f.sub('S_TIMESTP','OUT'),'w'){|fp|
    fp << File.read(f).sub(/Fields=\s*\d+/,'Fields= 1').sub(/2\).*?Element/m,'Element')
  } unless f['OUT']
}

