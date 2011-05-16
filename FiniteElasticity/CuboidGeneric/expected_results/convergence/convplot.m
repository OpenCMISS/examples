%%
a=[384	5.47E-008	5.45E-008	7.11E-006	7.11E-006
1272	1.53E-008	1.37E-008	4.00E-007	4.00E-007
8176	2.79E-007	2.79E-007	2.56E-008	2.55E-008
58080	0.49636E-09 0.14520E-09   0.16946E-08   0.17042E-08];

loglog(a(:,1),a(:,4),'.-',a(:,1),a(:,5),'x-');
xlabel('Degrees of freedom');
ylabel('Integrated^2 Error');
legend('x variable','y variable');