
clear
clc

Y = 1.4 ;
%{
output = fopen('M-PM.txt','w') ;

% For interpolation of Prandtl-Meyer function and Mach number
fprintf(output,'M\t\tv\n') ;
for i = 1 : 0.001 : 5
    v_x = PM(i,Y) ;
    fprintf(output,'%f\t%f\n',i,v_x) ;
end
%}

data = readtable('M-PM.txt') ;
M_i = data{:,1} ;
v_i = data{:,2} ;

%nl = 2 ; % Number of characteristic lines
nl = 20 ; % Number of characteristic lines
nlines = nl ;
Me = 2 ;
%Me = 2.95 ;
%Me = 2.4 ;
Ae_At = E(Me,Y) ;

ve = PM(Me,Y) ;
theta_max = ve/2 ;
d_theta = theta_max / nl ;
theta_t = [] ;
v_t = [] ;
for a = 1:1:nl
    theta_t(a) = d_theta * a ;
    v_t(a) = theta_t(a) ;
    M_t(a) = Mx(v_i,M_i,v_t(a)) ;
    u_t(a) = ux(M_t(a)) ;
end

x_t = 0 ;
%y_t = 1 ;
y_t = 1 ;
%y_t = 0.4811 ;

%{
% Fo testing
v_x = 26.3798 ;
M_x = Mx(v_i,M_i,v_x) ;
u_x = ux(M_x) ;
%}
%% Automation
% nl = nlines (or) number of lines
% np = npoints (or) number of points
% nc = ncentre (or) number of centreline points
% ni = ninner (or) number of inner points
% nw = nwall (or) number of wall points

np = nl + nl*(nl+1)/2 ;

j = 1 + nl ;
k = 0 ;

nw = [] ; % Wall points
nc = [] ; % Center points
nin = [] ; % Inner points

% Wall points
for i = 1 : 1 : np
    if i == j + k
        k = k + 1 ;
        nw(k) = i ;
        j = j + nl - k ;
    end
end

j = 1 ;
% Center points
while j < length(nw) 
nc(j) = nw(j) + 1 ;
j = j + 1 ;
end

nc = [1,nc] ;

% Inner points
point = 1 ;
for i = 1 : 1 : np
if ~ismember(i, nc) && ~ismember(i, nw)
nin(point) = i ;
point = point + 1 ;
end
end
a = 1 ; % a is starting from 2 as the value for the first throat point in already taken
b = 1 ;
c = 1 ;
d = 1 ;
e = 1 ;
f = 1 ;
x = [] ;
y = [] ;
theta = [] ;
v = [] ;
M = [] ;
u = [] ;
% Values for center points
x_c = [] ;
y_c = [] ;
theta_c = [] ;
v_c = [] ;
M_c = [] ;
m_minus_c = [] ;
tm_minus_c = [] ;
% Values for wall points
x_w = [] ;
y_w = [] ;
theta_w = [] ;
v_w = [] ;
M_w = [] ;
m_plus_w = [] ;
tm_plus_w = [] ;
m_plus_i_w = [] ;
tm_plus_i_w = [] ;
% values for inner points
x_in = [] ;
y_in = [] ;
theta_in = [] ;
v_in = [] ;
M_in = [] ;
u_in = [] ;
m_plus_in = [] ;
theta_in = [] ;
w = 0 ;
z = 1 + nl ;
syms xx yy
%{
% Re-intialization of j and k
j = 1 + nl ;
k = 0 ;
%}
for nl1 = 1:1:length(nw)
    for nl2 = nc(nl1):1:nw(nl1)
        for i = nl2
            %% CENTER POINTS
            if (a <= length(nc)) && i == nc(a)
                % 1ST POINT
                if i == 1
                    theta(i) = 0 ;
                    theta_c(a) = theta(i) ;
                    v(i) = K_minus(theta_t(i),v_t(i)) ;
                    M(i) = Mx(v_i,M_i,v(i)) ; % Corresponding Mach
                    u(i) = ux(M(i)) ; % Corresponding Mach angle
                    m_minus_t_1 = m_minus(theta_t(i),theta(i),u_t(i),u(i)) ; % Slope angle
                    tm_minus_t_1 = t(m_minus_t_1) ; % Slope
                    y(i) = 0 ;
                    %x(i) = -1 / tm_minus_t_1 ;
                    eq1 = tm_minus_t_1 * (xx - x_t) + y_t == 0 ;
                    s = solve(eq1,xx) ;
                    x(i) = s ;
                    x_c(a) = x(i) ;
                    y_c(a) = y(i) ;
            
                elseif i == nc(b+1) && (b+1 <= length(nc))
                    % CENTER POINTS
                    theta(i) = 0 ;
                    theta_c(b+1) = theta(i) ;
                    v(i) = K_minus(theta(nc(b)+1),v(nc(b)+1)) ;
                    M(i) = Mx(v_i,M_i,v(i)) ; % Corresponding Mach
                    u(i) = ux(M(i)) ; % Corresponding Mach angle
                    m_minus_c(b+1) = m_minus(theta(nc(b)+1),theta(i),u(nc(b)+1),u(i)) ; % Slope angle
                    tm_minus_c(b+1) = t(m_minus_c(b+1)) ;
                    y(i) = 0 ;
                    eq1 = tm_minus_c(b+1) * (xx - x(nc(b)+1)) + y(nc(b)+1) == 0 ;
                    s = solve(eq1,xx) ;
                    x(i) = s ;%double(s) ;
                    x_c(a) = x(i) ;
                    y_c(a) = y(i) ;
                    b = b + 1 ;
                end
                a = a + 1 ;
            end

            %% INNER POINTS
            if (c <= length(nin)) && (i == nin(c))
                if i < nw(1)
                    theta(i) = ( K_plus(theta(i-1),v(i-1)) + K_minus(theta_t(d+1),v_t(d+1)) ) / 2 ;
                    theta_in(d) = theta(i) ; 
                    v(i) = ( K_minus(theta_t(d+1),v_t(d+1)) - K_plus(theta(i-1),v(i-1)) ) / 2 ;
                    v_in(d) = v(i) ;
                    
                    M(i) = Mx(v_i,M_i,v(i)) ; % Corresponding Mach
                    M_in(d) = M(i) ;
                    u(i) = ux(M(i)) ; % Corresponding Mach angle
                    u_in(d) = u(i) ;
                    m_plus_in(d) = m_plus(theta(i-1),theta(i),u(i-1),u(i)) ;
                    tm_plus_in(d) = t(m_plus_in(d)) ;
                    m_minus_in(d) = m_minus(theta_t(d+1),theta(i),u_t(d+1),u(i)) ;
                    tm_minus_in(d) = t(m_minus_in(d)) ;
                    
                    eqn1 = yy - y(i-1) == tm_plus_in(d)*(xx-x(i-1)) ;
                    eqn2 = yy - y_t == tm_minus_in(d)*(xx-x_t) ;
                    s = solve([eqn1,eqn2],[xx,yy]) ;
                    x(i) = s.xx ;%double(s.xx) ;
                    y(i) = s.yy ;%double(s.yy) ;
                    
                    x_in(c) = x(i) ;
                    y_in(c) = y(i) ;
                    
                    d = d + 1 ;
                elseif i < nw(e)
                    theta(i) = ( K_plus(theta(i-1),v(i-1)) + K_minus(theta(i-nl),v(i-nl)) ) / 2 ;
                    theta_in(d) = theta(i) ; 
                    v(i) = ( K_minus(theta(i-nl),v(i-nl)) - K_plus(theta(i-1),v(i-1)) ) / 2 ;
                    v_in(d) = v(i) ;
                    
                    M(i) = Mx(v_i,M_i,v(i)) ; % Corresponding Mach
                    M_in(d) = M(i) ;
                    u(i) = ux(M(i)) ; % Corresponding Mach angle
                    u_in(d) = u(i) ;
                    m_plus_in(d) = m_plus(theta(i-1),theta(i),u(i-1),u(i)) ;
                    tm_plus_in(d) = t(m_plus_in(d)) ;
                    m_minus_in(d) = m_minus(theta_in(d),theta(i-nl),u_in(d),u(i-nl)) ;
                    tm_minus_in(d) = t(m_minus_in(d)) ;
                    
                    eqn1 = yy - y(i-1) == tm_plus_in(d)*(xx-x(i-1)) ;
                    eqn2 = yy - y_t == tm_minus_in(d)*(xx-x_t) ;
                    s = solve([eqn1,eqn2],[xx,yy]) ;
                    x(i) = s.xx ;%double(s.xx) ;
                    y(i) = s.yy ;%double(s.yy) ;
                    
                    x_in(c) = x(i) ;
                    y_in(c) = y(i) ;
                    
                    d = d + 1 ;
                end
                %{
                Is it possible to make a list in matlab for inner points ,
                => For 3 characteristic lines , i want a list looking something 
                like nin = [ (2,3) , (6) ]
                => For 3 characteristic lines , i want a list looking something 
                like nin = [ (2,3,4) , (7,8) , (11) ]
                %}
                
                c = c + 1 ;
            end
            %% WALL POINTS
            if i == nw(e)
                v(i) = v(i-1) ;
                M(i) = Mx(v_i,M_i,v(i)) ; % Corresponding Mach
                u(i) = ux(M(i)) ; % Corresponding Mach angle
            
                if i == nw(1)
                    theta(i) = theta(i-1) ;
                    theta_w(f) = theta(i) ;
                    m_plus_w(f) = m_plus(theta(i-1),theta(i),u(i-1),u(i)) ;
                    tm_plus_w(f) = t(m_plus_w(f)) ;
                    m_t(f) = (theta_max + theta(i)) / 2 ; % Slope angle from throat
                    tm_t_w(f) = t(m_t(f)) ;
                    
                    eqn1 = yy - y(i-1) == tm_plus_w(f)*(xx-x(i-1)) ;
                    eqn2 = yy - y_t == tm_t_w(f)*(xx-x_t) ;
                    s = solve([eqn1,eqn2],[xx,yy]) ;
                    x(i) = s.xx ;%double(s.xx) ;
                    y(i) = s.yy ;%double(s.yy) ;
                    x_w(f) = x(i) ;
                    y_w(f) = y(i) ;
                    
                    f = f + 1 ;
                elseif (i > nw(1)) && (i < nw(end))
                    theta(i) = theta(i-1) ;
                    theta_w(f) = theta(i) ;
                    m_plus_w(f) = (theta_w(f-1) + theta(i)) / 2 ;
                    tm_plus_w(f) = t(m_plus_w(f)) ;
                    m_plus_i_w(f) = m_plus(theta(i-1),theta(i),u(i-1),u(i)) ;
                    tm_plus_i_w(f) = t(m_plus_i_w(f)) ;
                    
                    eqn1 = yy - y(i-1) == tm_plus_i_w(f)*(xx-x(i-1)) ;
                    eqn2 = yy - y_w(f-1) == tm_plus_w(f)*(xx-x_w(f-1)) ;
                    s = solve([eqn1,eqn2],[xx,yy]) ;
                    x(i) = s.xx ;%double(s.xx) ;
                    y(i) = s.yy ;%double(s.yy) ;
                    x_w(f) = x(i) ;
                    y_w(f) = y(i) ;
                    
                    f = f + 1 ;
                elseif i == nw(end)
                    theta(i) = 0 ;
                    theta_w(f) = theta(i) ;
                    m_plus_w(f) = (theta_w(f-1) + theta(i)) / 2 ;
                    tm_plus_w(f) = t(m_plus_w(f)) ;
                    m_plus_i_w(f) = m_plus(theta(i-1),theta(i),u(i-1),u(i)) ;
                    tm_plus_i_w(f) = t(m_plus_i_w(f)) ;
                    
                    eqn1 = yy - y(i-1) == tm_plus_i_w(f)*(xx-x(i-1)) ;
                    eqn2 = yy - y_w(f-1) == tm_plus_w(f)*(xx-x_w(f-1)) ;
                    s = solve([eqn1,eqn2],[xx,yy]) ;
                    x(i) = s.xx ;%double(s.xx) ;
                    y(i) = s.yy ;%double(s.yy) ;
                    x_w(f) = x(i) ;
                    y_w(f) = y(i) ;
                    
                    f = f + 1 ;
                end
                e = e + 1 ;
            end
        end
    end
if i > nw(1)
    nl = nl - 1 ;
end
%nl1 = nw(nl1);
end
%% PLOTTING
clc
Ae_At_code = ( pi*y(end)^2 )/( pi*y_t^2 ) ;
De_Dt = ( 2*y(end) )/( 2*y_t ) ;
fprintf("Ae/At = %f\nDe/Dt = %f\n",Ae_At_code,De_Dt)

hold on 
% Axis of Symmetry
line([x_w(end),x_t],[0,0],Linestyle='--',Linewidth=1,color='black')

% Wall contour
line([x_t,x_w(1)],[y_t,y_w(1)],Linewidth=2,color='black') ;
for i = 2:1:length(x_w)
    line([x_w(i-1),x_w(i)],[y_w(i-1),y_w(i)],Linewidth=2,color='black') ;
end

%{
% Straight lines
for i = 1:1:length(x_w)
line([x_t,x_c(i)],[y_t,y_c(i)],Linewidth=2,color='red') ;
line([x_w(i),x_c(i)],[y_w(i),y_c(i)],Linewidth=2,color='red') ;
end
%}

% Points
for i = 1:1:np
    plot(x(i),y(i),marker='.',color='blue',MarkerSize=15) ;
end

%{
% Connecting lines
cl1 = 1 ;
cl2 = 1 ;
cl3 = 0 ;
cl4 = 1 ;
cl5 = 1 ;
for cl1 = 1:1:nlines
    for cl2 = nc(cl1):1:nw(cl1)
        for i = cl2
            % CENTERS
            if i == nc(cl1)
            % THROAT TO CENTER POINTS
                if i == 1
                line([x_t,x_c(1)],[y_t,y_c(1)],Linewidth=2,color='red') ;
                
                % INNER TO CENTER POINTS
                elseif i == nw(cl1) && i > 1
                %line([x_in(cl3),x_c(cl3)],[y_in(cl3),y_c(cl3)],Linewidth=2,color='red') ;
                line([x(x_c(cl3)+1),x_c(cl1)],[y(y_c(cl3)+1),y_c(cl1)],Linewidth=2,color='red') ;
                end
            end
            % INNERS
            if i == nin(cl4)
                % THROAT TO INNER
                if i < nw(1)
                line([x_t,x_in(cl4)],[y_t,y_in(cl4)],Linewidth=2,color='red') ;
                cl4 = cl4 + 1 ;
                
                % CENTER TO INNER
                elseif i == nc(cl1)+1
                line([x_c(cl1),x_in(cl5)],[y_c(cl1),y_in(cl5)],Linewidth=2,color='red') ;
                cl5 = cl5 + 1 ;
                % INNER TO INNER
                elseif i <= nw(cl1)
                %{
                if i == nin()
                cl4 = cl6 + 1 ;
                end
                %}
                line([x_in(cl6-1),x_in(cl6)],[y_in(cl6-1),y_in(cl6)],Linewidth=2,color='red') ;
                %{
                % INNER TO WALL
                elseif i == nw(cl1)
                line([x(i),x_in(cl4)],[y(i),y_in(cl4)],Linewidth=2,color='red') ;
                %}
                end
            end
        end
    end
    cl3 = cl3 + 1 ;
end
%}


cl1 = 1 ;
cl2 = 1 ;
cl3 = 1 ;
cl4 = -1 ;
for cl1 = 1:1:nlines
    for cl2 = nc(cl1):1:nw(cl1)

        % THROAT TO CENTER
        %line([x_t,x_c(cl1)],[y_t,y_c(cl1)],Linewidth=2,color='red') ;
        if cl1 == 1
            line([x_t,x_c(cl1)],[y_t,y_c(cl1)],Linewidth=2,color='red') ;
        elseif cl1 > 1
            line([x_t,x_in(cl1-1)],[y_t,y_in(cl1-1)],Linewidth=2,color='red') ;
        end
        
        % CENTER TO INNER
        if cl2 == nc(cl3) + 1
            line([x_c(cl1),x(nc(cl1)+1)],[y_c(cl1),y(nc(cl1)+1)],Linewidth=2,color='red') ;
            cl3 = cl3 + 1 ;
        elseif cl2 > nc(cl1) && cl2 <= nw(cl1)
            line([x(cl2),x(cl2-1)],[y(cl2),y(cl2-1)],Linewidth=2,color='red') ;
        end

        % INNER TO INNER
        if cl1 > 1 && cl2 < nw(cl1)
            line([x(cl2-(nlines-cl4)),x(cl2)],[y(cl2-(nlines-cl4)),y(cl2)],Linewidth=2,color='red') ;
        end
    end
    cl4 = cl4 + 1 ;
end

axis equal

%% FUNCTIONS

function v = PM(M,Y)
    v = sqrt( (Y+1)/(Y-1) ) * atand( sqrt( ((Y-1)/(Y+1))*(M^2 - 1) ) ) - atand( sqrt(M^2 - 1) ) ;
end

% Mach number and Prandtl-Meyer function
function S = Mx(v1,M1,v2)
    S = interp1(v1,M1,v2) ;
end

% Mach angle
function u = ux(M)
    u = asind(1/M) ;
end

% K- characteristic line
function z = K_minus(theta1,v1)
    z = theta1 + v1 ;
end

% K+ characteristic line
function z = K_plus(theta1,v1)
    z = theta1 - v1 ;
end

% Angle of C- characteristic slope
function s = m_minus(Ox,Oy,ux,uy)
    s = ( (Ox-ux) + (Oy-uy) ) / 2 ;
end

% Angle of C+ characteristic slope
function s = m_plus(Ox,Oy,ux,uy)
    s = ( (Ox+ux) + (Oy+uy) ) / 2 ;
end

% Slope
function s = t(u)
    s = tand(u) ;
end

% Area ratio
function s = E(M,Y)
    a = (Y+1)/2 ;
    b = (Y-1)/2 ;
    s = (1/M) * ((1/a) * (1 + (b*(M^2))))^(a/(2*b));
end