%% scaling


% Generated on a single proc, 8-core shared memory machine
% using direct solver (MUMPS)
p=[1 2 4 8];
t1=[426.474 250.853 146.534 92.884];
t2=[13*60+2.66 7*60+24.76 4*60+5.886 2*60+19.3];
t3=[26*60+32.3 14*60+44.899 8*60+30.869 4*60+44.726];
ratio1=t1(1)./t1;
ratio2=t2(1)./t2;
ratio3=t3(1)./t3;
plot(p,ratio1,'.-',p,ratio2,'rx-',p,ratio3,'kd-');
legend('DOF=6144','DOF=12288','DOF=19968','location','northwest');
xlabel('number of processors');
ylabel('speedup');