s1 = sin(2*pi/5);
s2 = sin(4*pi/5);
c1 = cos(2*pi/5);
c2 = cos(pi/5);
V = [0 1;-s1 c1;-s2 -c2;s2 -c2;s1 c1];
d12 = sqrt(2*(1-c1));
% d12 = d15 = d23 = d34 = d45
d13 = sqrt(2*(1+c2));
% d13 = d14
delta = 1;
kv = 10;
T = 0.01;

q0 = zeros(5,2);
qout = zeros(1000,5,2);
for i = 1:5
    q0(i,:) = [V(i,:)+delta*(rand(1,2)-0.5*ones(1,2))];
end
q = q0;

% plot(V(:,1),V(:,2),'o')
% hold on
% plot(q0(:,1),q0(:,2),'*')

for t = 1:1000
    pt = 1+0.5*sin(0.4*t*T);
    pt_dot = 0.2*cos(0.4*t*T);
    dij = pt*[d12 d13 d13 d12 d12 d12 d12];
    dij_dot = pt_dot*[d12 d13 d13 d12 d12 d12 d12];
    dv = dij.*dij_dot;
    
    vd = [1 cos(t*T)];

    Rq = [(q(1,:)-q(2,:))' (q(2,:)-q(1,:))' zeros(2,3);
          (q(1,:)-q(3,:))' zeros(2,1) (q(3,:)-q(1,:))' zeros(2,2);
          (q(1,:)-q(4,:))' zeros(2,2) (q(4,:)-q(1,:))' zeros(2,1);
          (q(1,:)-q(5,:))' zeros(2,3) (q(5,:)-q(1,:))';
          zeros(2,1) (q(2,:)-q(3,:))' (q(3,:)-q(2,:))' zeros(2,2);
          zeros(2,2) (q(3,:)-q(4,:))' (q(4,:)-q(3,:))' zeros(2,1);
          zeros(2,3) (q(4,:)-q(5,:))' (q(5,:)-q(4,:))'];
    
    Rq_plus = Rq'/(Rq*Rq');

    qij = [norm((q(1,:)-q(2,:))) norm((q(1,:)-q(3,:))) norm((q(1,:)-q(4,:))) ...
        norm((q(1,:)-q(5,:))) norm((q(2,:)-q(3,:))) norm((q(3,:)-q(4,:))) norm((q(4,:)-q(5,:)))];

    e = qij - dij;

    z = e.*(e+2*dij);

    u = zeros(5,2);
    Rqt = Rq_plus;
    for i = 1:5
        k = (-kv*z+dv);
        Rqi = Rqt(i,1:2)*k(1) + Rqt(i,3:4)*k(2) + Rqt(i,5:6)*k(3) + Rqt(i,7:8)*k(4) ...
            + Rqt(i,9:10)*k(5) + Rqt(i,11:12)*k(6) + Rqt(i,13:14)*k(7);
        u(i,:) = Rqi + vd;
    end

    q = q+u*T;
    qout(t,:,:) = q;
end

figure(1)
plot(qout(1,:,1),qout(1,:,2),'.');
hold on
grid on
for t = 2:1000
    plot(qout(t,:,1),qout(t,:,2),'.');
end
plot(qout(end,:,1),qout(end,:,2),'o');
plot(q0(:,1),q0(:,2),'*');
hold off
