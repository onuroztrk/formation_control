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
e_out = zeros(1000,9);
qij_out = zeros(1000,9);
for i = 1:5
    q0(i,:) = [V(i,:)+delta*(rand(1,2)-0.5*ones(1,2))];
end
q = q0;

for t = 1:1000

    Rq = [(q(1,:)-q(2,:)) (q(2,:)-q(1,:)) zeros(1,8);
          (q(1,:)-q(6,:)) zeros(1,8) (q(6,:)-q(1,:));
          zeros(1,2) (q(2,:)-q(3,:)) (q(3,:)-q(2,:)) zeros(1,6);
          zeros(1,2) (q(2,:)-q(6,:)) zeros(1,6) (q(6,:)-q(2,:));
          zeros(1,4) (q(3,:)-q(4,:)) (q(4,:)-q(3,:)) zeros(1,4);
          zeros(1,4) (q(3,:)-q(6,:)) zeros(1,4) (q(6,:)-q(3,:));
          zeros(1,6) (q(4,:)-q(5,:)) (q(5,:)-q(4,:)) zeros(1,2);
          zeros(1,6) (q(4,:)-q(6,:)) zeros(1,2) (q(6,:)-q(4,:));
          zeros(1,8) (q(5,:)-q(6,:)) (q(6,:)-q(5,:))];

    qij = [norm((q(1,:)-q(2,:))) norm((q(1,:)-q(6,:))) norm((q(2,:)-q(3,:))) norm((q(2,:)-q(6,:))) norm((q(3,:)-q(4,:))) ...
            norm((q(3,:)-q(6,:))) norm((q(4,:)-q(5,:))) norm((q(4,:)-q(6,:))) norm((q(5,:)-q(6,:)))];
        
    qij_out(t,:) = qij;
    
    qin = [(q(1,:)-q(6,:));(q(2,:)-q(6,:));(q(3,:)-q(6,:));(q(4,:)-q(6,:));(q(5,:)-q(6,:));0 0];
    e = qij - dij;
    e_out(t,:) = e;

    z = e.*(e+2*dij);
    v0 = [1 cos(t*T)];
    u = zeros(6,2);
    Rqt = Rq';
    for i = 1:6
        v = cross(w0,[qin(i,:) 0]);
        sum = 0;
        for j = 1:9
            sum = sum + Rqt(2*(i-1)+1:2*(i-1)+2,j)*z(j);
        end
        u(i,:) = -kv*sum'+v0+v(1:2);
    end

    q = q+u*T;
    qout(t,:,:) = q;
end

g1 = graph([1 1 2 2 3 3 4 4 5],[2 6 3 6 4 6 5 6 6],qij_out(1,:));
plot(g1,'XData',qout(1,:,1)','YData',qout(1,:,2)');
hold on
grid on

% for t = 2:1000
%     plot(qout(t,:,1),qout(t,:,2),'.');
% end
for t = 200:200:1000
    g1 = graph([1 1 2 2 3 3 4 4 5],[2 6 3 6 4 6 5 6 6],qij_out(t,:));
    plot(g1,'XData',qout(t,:,1)','YData',qout(t,:,2)');
end
hold off

figure(2)
plot(0:0.01:0.99,e_out(1:100,:))
ylabel('e_{ij}')
xlabel('s')
grid on