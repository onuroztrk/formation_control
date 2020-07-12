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
qout = zeros(1500,5,2);
e_out = zeros(1500,7);
for i = 1:5
    q0(i,:) = [V(i,:)+delta*(rand(1,2)-0.5*ones(1,2))];
end
q = q0;

for t = 1:1500
    pt = 1+0.5*sin(0.4*t*T);
    dij = pt*[d12 d13 d13 d12 d12 d12 d12];

    if t == 1
        dij_ = dij;
    end
    ddij = (dij-dij_)/T;
    dij_ = dij;
    dv = dij.*ddij;
    vd = [1 cos(t*T)];

    Rq = [(q(1,:)-q(2,:)) (q(2,:)-q(1,:)) zeros(1,6);
          (q(1,:)-q(3,:)) zeros(1,2) (q(3,:)-q(1,:)) zeros(1,4);
          (q(1,:)-q(4,:)) zeros(1,4) (q(4,:)-q(1,:)) zeros(1,2);
          (q(1,:)-q(5,:)) zeros(1,6) (q(5,:)-q(1,:));
          zeros(1,2) (q(2,:)-q(3,:)) (q(3,:)-q(2,:)) zeros(1,4);
          zeros(1,4) (q(3,:)-q(4,:)) (q(4,:)-q(3,:)) zeros(1,2);
          zeros(1,6) (q(4,:)-q(5,:)) (q(5,:)-q(4,:))];
    
    Rq_plus = Rq'/(Rq*Rq');

    qij = [norm((q(1,:)-q(2,:))) norm((q(1,:)-q(3,:))) norm((q(1,:)-q(4,:))) ...
        norm((q(1,:)-q(5,:))) norm((q(2,:)-q(3,:))) norm((q(3,:)-q(4,:))) norm((q(4,:)-q(5,:)))];

    e = qij - dij;
    e_out(t,:) = e;
    z = e.*(e+2*dij);

    u = zeros(5,2);
    Rqt = Rq_plus;
    for i = 1:5
        k = (-kv*z+dv);
        sum = 0;
        for j = 1:7
            sum = sum + Rqt(2*(i-1)+1:2*(i-1)+2,j)*k(j);
        end
        u(i,:) = sum' + vd;
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
plot(qout(1000,:,1),qout(800,:,2),'o');
plot(qout(1200,:,1),qout(800,:,2),'o');
plot(qout(1400,:,1),qout(800,:,2),'o');
plot(qout(end,:,1),qout(end,:,2),'o');
plot(q0(:,1),q0(:,2),'*');
hold off

figure(2)
plot(0:0.01:0.99,e_out(1:100,:))
ylabel('e_{ij}')
xlabel('s')
grid on
