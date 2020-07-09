s1 = sin(2*pi/5);
s2 = sin(4*pi/5);
c1 = cos(2*pi/5);
c2 = cos(pi/5);
V = [0 1;-s1 c1;-s2 -c2;s2 -c2;s1 c1];
d12 = sqrt(2*(1-c1));
% d12 = d15 = d23 = d34 = d45
d13 = sqrt(2*(1+c2));
% d13 = d14
dij = [d12 d13 d13 d12 d12 d12 d12];
delta = 1;
kv = 1;
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

    Rq = [(q(1,:)-q(2,:)) (q(2,:)-q(1,:)) zeros(1,6);
          (q(1,:)-q(3,:)) zeros(1,2) (q(3,:)-q(1,:)) zeros(1,4);
          (q(1,:)-q(4,:)) zeros(1,4) (q(4,:)-q(1,:)) zeros(1,2);
          (q(1,:)-q(5,:)) zeros(1,6) (q(5,:)-q(1,:));
          zeros(1,2) (q(2,:)-q(3,:)) (q(3,:)-q(2,:)) zeros(1,4);
          zeros(1,4) (q(3,:)-q(4,:)) (q(4,:)-q(3,:)) zeros(1,2);
          zeros(1,6) (q(4,:)-q(5,:)) (q(5,:)-q(4,:))];

    qij = [norm((q(1,:)-q(2,:))) norm((q(1,:)-q(3,:))) norm((q(1,:)-q(4,:))) ...
        norm((q(1,:)-q(5,:))) norm((q(2,:)-q(3,:))) norm((q(3,:)-q(4,:))) norm((q(4,:)-q(5,:)))];

    e = qij - dij;

    z = e.*(e+2*dij);

    u = zeros(5,2);
    Rqt = Rq';
    for i = 1:5
        sum = 0;
        for j = 1:7
            sum = sum + Rqt(2*(i-1)+1:2*(i-1)+2,j)*z(j);
        end
        u(i,:) = -kv*sum';
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
