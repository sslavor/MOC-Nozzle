output = fopen('M-PM.txt','w') ;
Y = 1.4 ;

% For interpolation of Prandtl-Meyer function and Mach number
fprintf(output,'M\t\tv\n') ;

for i = 1 : 1e-5 : 5
    v_x = PM(i,Y) ;
    fprintf(output,'%f\t%f\n',i,v_x) ;
end

function v = PM(M,Y)
    v = sqrt( (Y+1)/(Y-1) ) * atand( sqrt( ((Y-1)/(Y+1))*(M^2 - 1) ) ) - atand( sqrt(M^2 - 1) ) ;
end