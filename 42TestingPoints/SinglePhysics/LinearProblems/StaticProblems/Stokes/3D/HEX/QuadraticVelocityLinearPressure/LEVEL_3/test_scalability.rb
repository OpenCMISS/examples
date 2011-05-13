mpi_path = '${OPENCMISSEXTRAS_ROOT}/cm/external/x86_64-linux/mpich2/gnu/bin'
system("#{mpi_path}/mpd&")

ns = [1,2,4,8]
summary = {}

File.readlines('tests').each{|line|
  cmd = line.strip
  next if cmd.empty? || cmd[0..0]=='#'
  cmd = cmd.split(' ')
  Dir.chdir(cmd[0]){
  `make allclean; make opt64 COMPILER=gnu` # compile opt.
  ns.each{|n|
  `rm -f time.tmp`
  fullcmd = "./bin/x86_64-linux/mpich2/gnu/#{cmd[0]}Example #{cmd[1..-1].join(' ')}"
  mpicmd = "#{mpi_path}/mpiexec -n #{n} /usr/bin/time -a -o time.tmp -f \"%e %U %S %M %x\" #{fullcmd} &> /dev/null"
  `#{mpicmd}`
  puts '='*50, mpicmd, "#{n} PROCESSORS", '-'*20
  data = File.readlines('time.tmp').map{|l| 
    w,u,s,m,x = l.split(' ')
      puts "Walltime #{w} User #{u} System #{s} Memory #{m}  Exit code #{x}"
    [w,u+s,m,x]
  }
  summary[fullcmd] ||= []
  summary[fullcmd] << [n] + data.transpose.map{|d| d.inject{|a,b|a.to_f+b.to_f}.to_f/d.size }
  summary[fullcmd][-1][-1] += 1 unless data.size==n # force error when unexpected number of processors
 }
 }
}
 puts "=============================== SUMMARY =============================== "

 File.open('summary.out','w'){|fo|
  [fo,$stdout].each{|f|
   summary.sort.each{|c,d|
    f.puts '-'*50,"COMMAND: #{c}", "Processors\tWalltime\tCPU Time\tMemory"
    d.each{|data| n,w,t,m,x = *data; f.puts "#{n}       \t#{'%.2f'%(w/60)}m     \t#{'%.2f'%(t/60)}m     \t#{'%.2f'%(m/1048576)}Gb/cpu  \t#{x.to_i.zero? ? 'OK' : 'ERROR!'}" } rescue puts "<ERROR PROCESSING DATA> #{d.inspect}"
   }
  }
 }




