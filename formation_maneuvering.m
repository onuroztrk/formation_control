w0 = [0 0 1];
s1 = sin(2*pi/5);
s2 = sin(4*pi/5);
c1 = cos(2*pi/5);
c2 = cos(pi/5);
V = [0 1;-s1 c1;-s2 -c2;s2 -c2;s1 c1; 0 0];
d12 = sqrt(2*(1-c1));
% d12 = d23 = d34 = d45
d16 = 1;
% d16 = d26 = d36 = d46 = d56
dij = [d12 d16 d12 d16 d12 d16 d12 d16 d16];
delta = 1;
kv = 1;
T = 0.01;

q0 = zeros(6,2);
qout = zeros(1000,6,2);
for i = 1:5
    q0(i,:) = [V(i,:)+delta*(rand(1,2)-0.5*ones(1,2))];
end
q = q0;

for t = 1:1000

    Rq = [(q(1,:)-q(2,:))' (q(2,:)-q(1,:))' zeros(2,4);
          (q(1,:)-q(6,:))' zeros(2,4) (q(6,:)-q(1,:))';
          zeros(2,1) (q(2,:)-q(3,:))' (q(3,:)-q(2,:))' zeros(2,3);
          zeros(2,1) (q(2,:)-q(6,:))' zeros(2,3) (q(6,:)-q(2,:))';
          zeros(2,2) (q(3,:)-q(4,:))' (q(4,:)-q(3,:))' zeros(2,2);
          zeros(2,2) (q(3,:)-q(6,:))' zeros(2,2) (q(6,:)-q(3,:))';
          zeros(2,3) (q(4,:)-q(5,:))' (q(5,:)-q(4,:))' zeros(2,1);
          zeros(2,3) (q(4,:)-q(6,:))' zeros(2,1) (q(6,:)-q(4,:))';
          zeros(2,4) (q(5,:)-q(6,:))' (q(6,:)-q(5,:))'];

    qij = [norm((q(1,:)-q(2,:))) norm((q(1,:)-q(6,:))) norm((q(2,:)-q(3,:))) norm((q(2,:)-q(6,:))) norm((q(3,:)-q(4,:))) ...
            norm((q(3,:)-q(6,:))) norm((q(4,:)-q(5,:))) norm((q(4,:)-q(6,:))) norm((q(5,:)-q(6,:)))];
    
    qin = [(q(1,:)-q(6,:));(q(2,:)-q(6,:));(q(3,:)-q(6,:));(q(4,:)-q(6,:));(q(5,:)-q(6,:));0 0];
    e = qij - dij;

    z = e.*(e+2*dij);
    v0 = [1 cos(t*T)];
    u = zeros(6,2);
    Rqt = Rq';
    for i = 1:6
        v = cross(w0,[qin(i,:) 0]);
        u(i,:) = -kv*(Rqt(i,1:2)*z(1)+Rqt(i,3:4)*z(2)+Rqt(i,5:6)*z(3)+Rqt(i,7:8)*z(4)+Rqt(i,9:10)*z(5)+ ...
            Rqt(i,11:12)*z(6)+Rqt(i,13:14)*z(7)+Rqt(i,15:16)*z(8)+Rqt(i,17:18)*z(9))+v0+v(1:2);
    end

    q = q+u*T;
    qout(t,:,:) = q;
end

figure(1)
plot(qout(1,:,1),qout(1,:,2),'.');
hold on
grid on
% for t = 2:1000
%     plot(qout(t,:,1),qout(t,:,2),'.');
% end
plot(qout(200,:,1),qout(200,:,2),'o');
plot(qout(400,:,1),qout(400,:,2),'o');
plot(qout(600,:,1),qout(600,:,2),'o');
plot(qout(800,:,1),qout(800,:,2),'o');
plot(qout(end,:,1),qout(end,:,2),'o');
plot(q0(:,1),q0(:,2),'*');
hold off