mpi_path = '${OPENCMISSEXTRAS_ROOT}/cm/external/x86_64-linux/mpich2/gnu/bin'
system("#{mpi_path}/mpd&")

ns = [1,2,4,8]
summary = {}

File.readlines('tests').each{|line|
  cmd = line.strip
  next if cmd.empty?
  cmd = cmd.split(' ')
  Dir.chdir(cmd[0]){
  `make allclean; make opt64 COMPILER=gnu` # compile opt.
  ns.each{|n|
  `rm -f time.tmp`
  fullcmd = "./bin/x86_64-linux/mpich2/gnu/#{cmd[0]}Example #{cmd[1..-1].join(' ')}"
  `#{mpi_path}/mpiexec -n #{n} /usr/bin/time -a -o time.tmp -f "%e %U %S %M %x" #{fullcmd} &> /dev/null`
  puts '-'*50, fullcmd, "#{n} PROCESSORS", '-'*50
  data = File.readlines('time.tmp').map{|l| 
    w,u,s,m,x = l.split(' ')
      puts "Walltime #{w} User #{u} System #{s} Memory #{m}  Exit code #{x}"
    [w,u+s,m,x]
  }
  summary[fullcmd] ||= []
  summary[fullcmd] << [n] + data.transpose.map{|d| d.inject{|a,b|a.to_f+b.to_f}.to_f/n }
 }
 }
}
 puts "------------------ SUMMARY ----------------"
  summary.each{|c,d|
    puts "COMMAND: #{c}", "Processors\tWalltime\tCPU Time\tMemory"
    d.each{|data| n,w,t,m,x = data.map{|s| s.to_s.ljust 10 }; puts "#{n}\t#{w}\t#{t}\t#{m}\t#{x.to_i.zero? ? 'OK' : 'ERRROR!'}" }
  }
 File.open('summary.out','w'){|f|
  summary.each{|c,d|
    f.puts "COMMAND: #{c}", "Processors\tWalltime\tCPU Time\tMemory"
    d.each{|data| n,w,t,m,x = data.map{|s| s.to_s.ljust 10 }; f.puts "#{n}\t#{w}\t#{t}\t#{m}\t#{x.to_i.zero? ? 'OK' : 'ERRROR!'}" }
  }
 }




