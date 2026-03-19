% First Order Control Barrier Function.
% The system has the form
% \dot{x} = u.
% The barrier function is defined as the norm between the system position
% and the obstacle position minus a raduis of 0.5.
% \norm(p_system - p_obstacle) - radius.

%% System Ax + Bu
I     = eye(2);
A     = [zeros(size(I))];
B     = I; 
K     = 1; %Controller gain
dt    = 0.001;  %Time step
t     = 20;     %Sim time
s     = 0:dt:t; %Time
xp = zeros(size(A,1),length(s));
xp(:,1) = [0; 0];    %Initial condotions
e = zeros(size(A,1),length(s)); %Error
r = zeros(size(A,1),length(s)); %Reference
%% Solver and control
for i = 1 : length(s)
    %Reference and Errors
    % r(:,i) = [3 + sin(0.5*(i-1)*dt); 5 + cos(0.5*(i-1)*dt)];
    r(:,i) = [3; 5];
    e(:,i) = xp(:,i) - r(:,i);

    %Barrier Function
    alfa1=5;
    hx1=norm([xp(1,i);xp(2,i)]-[1;2])-0.5;
    hx2=norm([xp(1,i);xp(2,i)]-[2.5;3])-0.5;
    gamma=eye(2);
    psi=[0;0];
    b1=([xp(1,i);xp(2,i)]-[1;2])'*gamma/norm([xp(1,i);xp(2,i)]-[1;2]);
    b2=([xp(1,i);xp(2,i)]-[2.5;3])'*gamma/norm([xp(1,i);xp(2,i)]-[2.5;3]);
    phi1=([xp(1,i);xp(2,i)]-[1;2])'*psi/norm([xp(1,i);xp(2,i)]-[1;2])+alfa1*hx1;
    phi2=([xp(1,i);xp(2,i)]-[2.5;3])'*psi/norm([xp(1,i);xp(2,i)]-[2.5;3])+alfa1*hx2;
    
    %Changing Barrier Function for a different obstacle
    hx=min([hx1,hx2]);
    b=[0 0]; phi=0;
    if hx==hx1
        b=b1;
        phi=phi1;
    elseif hx==hx2
        b=b2;
        phi=phi2;
    end

    %Nominal Controller
    ss1 = e(1,i);
    ss2 = e(2,i); 
    ux  = -K*(ss1);   %State feedback
    uy  = -K*(ss2);
    u_nominal = [ux; uy];
    
    %Barrier Activation Condition
    cond=phi+b*u_nominal;
    if cond<0
        u_safe=-b'*cond/norm(b)^2;
    else
        u_safe=[0;0];
    end
    
    %Control Barrier Function
    u = u_nominal + u_safe;
    
    %Numerical Integration
    xp(:,i+1) = xp(:,i) + dt*(A*xp(:,i) + B*u); %Forward Euler
end
%% Plots
%Errors (e1x, e2y)
figure(5)
grid on; hold on;
plot(s,e, 'LineWidth', 2);
xlabel({'Time'}, 'interpreter', 'latex', 'fontsize', 14)
ylabel({'Errors'}, 'interpreter', 'latex', 'fontsize', 14)
legend({'$e_{1}(t)$','$e_{2}(t)$'}, 'interpreter', 'latex', 'location', 'best', 'fontsize', 14);

%System Trajectory
figure(6)
% plot(3 + sin(0:pi/50:2*pi), 5 + cos(0:pi/50:2*pi))  %Trajectory
% hold on;
plot(1 + 0.5*sin(0:pi/50:2*pi), 2 + 0.5*cos(0:pi/50:2*pi))  %Obstacles
fill(1 + 0.5*sin(0:pi/50:2*pi), 2 + 0.5*cos(0:pi/50:2*pi), 'k')
hold on;
plot(2.5 + 0.5*sin(0:pi/50:2*pi), 3 + 0.5*cos(0:pi/50:2*pi))  %...
fill(2.5 + 0.5*sin(0:pi/50:2*pi), 3 + 0.5*cos(0:pi/50:2*pi), 'k')
grid on; hold on;
plot(xp(1,:),xp(2,:), 'LineWidth', 2)
axis([-1 5 -1 6]);
xlabel({'x'}, 'interpreter', 'latex', 'fontsize', 14)
ylabel({'y'}, 'interpreter', 'latex', 'fontsize', 14)
