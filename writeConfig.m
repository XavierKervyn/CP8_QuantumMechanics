%va écrire le configuration.in
filename = 'configuration.in';
fid = fopen(filename,'wt'); 

%Paramètres physiques:
fprintf(fid,['tfin=', num2str(tfin),'\n']);
fprintf(fid,['hbar=', num2str(hbar),'\n']);
fprintf(fid,['mass=', num2str(mass),'\n']); 
fprintf(fid,['xL=', num2str(xL),'\n']); 
fprintf(fid,['xR=', num2str(xR),'\n']); 
fprintf(fid,['omega=', num2str(omega),'\n']);
fprintf(fid,['Delta=', num2str(Delta_loc),'\n']);
fprintf(fid,['x0=', num2str(x0_loc),'\n']);
fprintf(fid,['n=', num2str(n),'\n']);
fprintf(fid,['sigma_norm=', num2str(sigma_norm),'\n']);

% Détecteur
fprintf(fid,['t_detect=', num2str(t_d_loc),'\n']);
fprintf(fid,['xda=', num2str(xda),'\n']);
fprintf(fid,['xdb=', num2str(xdb),'\n']);

%Paramètres numériques:
fprintf(fid,['dt=', num2str(dt_loc),'\n']);
fprintf(fid,['Ninters=', num2str(Ninters),'\n']);

%écriture dans des fichiers:
fprintf(fid,['output_potential=',fname2,'_potential.out\n']); 
fprintf(fid,['output_squared_wave = ',fname2,'_psi2.out\n']); 
fprintf(fid,['output_observables= ',fname2,'_observables.out\n']); 

fclose(fid);