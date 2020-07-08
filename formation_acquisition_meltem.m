V = [-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1;-1 -1 1;1 -1 1;1 1 1;-1 1 1];
d12 = 2;
d13 = 2*sqrt(2);
d14 = 2;
d15 = 2;
d18 = 2*sqrt(2);
d23 = 2;
d25 = 2*sqrt(2);
d26 = 2;
d34 = 2;
d36 = 2*sqrt(2);
d37 = 2;
d47 = 2*sqrt(2);
d48 = 2;
d56 = 2;
d58 = 2;
d67 = 2;
d68 = 2*sqrt(2);
d78 = 2;
dij = [d12 d13 d14 d15 d18 d23 d25 d26 d34 d36 d37 d47 d48 d56 d58 d67 d68 d78];

kv = 1;
ka = 1;
T = 0.01;

q0 = zeros(8,3);
v0 = zeros(8,3);
qout = zeros(1000,8,3);
for i = 1:8
    q0(i,:) = V(i,:)+(rand(1,3)-0.5*ones(1,3));
    v0(i,:) = 2*(rand(1,3)-0.5*ones(1,3));
end
q = q0;
q_ = q0;
vi = v0;
Rq_ = [(q(1,:)-q(2,:))' (q(2,:)-q(1,:))' zeros(3,6);
          (q(1,:)-q(3,:))' zeros(3,1) (q(3,:)-q(1,:))' zeros(3,5);
          (q(1,:)-q(4,:))' zeros(3,2) (q(4,:)-q(1,:))' zeros(3,4);
          (q(1,:)-q(5,:))' zeros(3,3) (q(5,:)-q(1,:))' zeros(3,3);
          (q(1,:)-q(8,:))' zeros(3,6) (q(8,:)-q(1,:))';
          zeros(3,1) (q(2,:)-q(3,:))' (q(3,:)-q(2,:))' zeros(3,5);
          zeros(3,1) (q(2,:)-q(5,:))' zeros(3,2) (q(5,:)-q(2,:))' zeros(3,3);
          zeros(3,1) (q(2,:)-q(6,:))' zeros(3,3) (q(6,:)-q(2,:))' zeros(3,2);
          zeros(3,2) (q(3,:)-q(4,:))' (q(4,:)-q(3,:))' zeros(3,4);
          zeros(3,2) (q(3,:)-q(6,:))' zeros(3,2) (q(6,:)-q(3,:))' zeros(3,2);
          zeros(3,2) (q(3,:)-q(7,:))' zeros(3,3) (q(7,:)-q(3,:))' zeros(3,1);
          zeros(3,3) (q(4,:)-q(7,:))' zeros(3,2) (q(7,:)-q(4,:))' zeros(3,1);
          zeros(3,3) (q(4,:)-q(8,:))' zeros(3,3) (q(8,:)-q(4,:))';
          zeros(3,4) (q(5,:)-q(6,:))' (q(6,:)-q(5,:))' zeros(3,2);
          zeros(3,4) (q(5,:)-q(8,:))' zeros(3,2) (q(8,:)-q(5,:))';
          zeros(3,5) (q(6,:)-q(7,:))' (q(7,:)-q(6,:))' zeros(3,1);
          zeros(3,5) (q(6,:)-q(8,:))' zeros(3,1) (q(8,:)-q(6,:))';
          zeros(3,6) (q(7,:)-q(8,:))' (q(8,:)-q(7,:))';];
      
for t = 1:1000

    Rq = [(q(1,:)-q(2,:))' (q(2,:)-q(1,:))' zeros(3,6);
          (q(1,:)-q(3,:))' zeros(3,1) (q(3,:)-q(1,:))' zeros(3,5);
          (q(1,:)-q(4,:))' zeros(3,2) (q(4,:)-q(1,:))' zeros(3,4);
          (q(1,:)-q(5,:))' zeros(3,3) (q(5,:)-q(1,:))' zeros(3,3);
          (q(1,:)-q(8,:))' zeros(3,6) (q(8,:)-q(1,:))';
          zeros(3,1) (q(2,:)-q(3,:))' (q(3,:)-q(2,:))' zeros(3,5);
          zeros(3,1) (q(2,:)-q(5,:))' zeros(3,2) (q(5,:)-q(2,:))' zeros(3,3);
          zeros(3,1) (q(2,:)-q(6,:))' zeros(3,3) (q(6,:)-q(2,:))' zeros(3,2);
          zeros(3,2) (q(3,:)-q(4,:))' (q(4,:)-q(3,:))' zeros(3,4);
          zeros(3,2) (q(3,:)-q(6,:))' zeros(3,2) (q(6,:)-q(3,:))' zeros(3,2);
          zeros(3,2) (q(3,:)-q(7,:))' zeros(3,3) (q(7,:)-q(3,:))' zeros(3,1);
          zeros(3,3) (q(4,:)-q(7,:))' zeros(3,2) (q(7,:)-q(4,:))' zeros(3,1);
          zeros(3,3) (q(4,:)-q(8,:))' zeros(3,3) (q(8,:)-q(4,:))';
          zeros(3,4) (q(5,:)-q(6,:))' (q(6,:)-q(5,:))' zeros(3,2);
          zeros(3,4) (q(5,:)-q(8,:))' zeros(3,2) (q(8,:)-q(5,:))';
          zeros(3,5) (q(6,:)-q(7,:))' (q(7,:)-q(6,:))' zeros(3,1);
          zeros(3,5) (q(6,:)-q(8,:))' zeros(3,1) (q(8,:)-q(6,:))';
          zeros(3,6) (q(7,:)-q(8,:))' (q(8,:)-q(7,:))';];
    
    vij = ((Rq-Rq_)/T)';
    Rq_ = Rq;

    qij = [norm((q(1,:)-q(2,:))) norm((q(1,:)-q(3,:))) norm((q(1,:)-q(4,:))) ...
           norm((q(1,:)-q(5,:))) norm((q(1,:)-q(8,:))) norm((q(2,:)-q(3,:))) ...
           norm((q(2,:)-q(5,:))) norm((q(2,:)-q(6,:))) norm((q(3,:)-q(4,:))) ...
           norm((q(3,:)-q(6,:))) norm((q(3,:)-q(7,:))) norm((q(4,:)-q(7,:))) ...
           norm((q(4,:)-q(8,:))) norm((q(5,:)-q(6,:))) norm((q(5,:)-q(8,:))) ...
           norm((q(6,:)-q(7,:))) norm((q(6,:)-q(8,:))) norm((q(7,:)-q(8,:)))];
       

    e = qij - dij;

    z = e.*(e+2*dij);

    u = zeros(8,3);
    Rqt = Rq';
    for i = 1:8
        sum = 0;
        for j = 1:3:54
             sum = sum + (ka*kv + 1)*Rqt(i,j:j+2)*z((j+2)/3) + (kv*(z((j+2)/3) + 2*Rqt(i,j:j+2)*Rqt(i,j:j+2)')*vij(i,j:j+2));
        end
        u(i,:) = -ka*vi(i,:) - sum;
    end

    vi = vi+u*T;
    q = q+vi*T;
    qout(t,:,:) = q;
end

figure(1)
plot3(qout(1,:,1),qout(1,:,2),qout(1,:,3),'.');
hold on
grid on
for t = 2:1000
    plot3(qout(t,:,1),qout(t,:,2),qout(t,:,3),'.');
end
plot3(qout(end,:,1),qout(end,:,2),qout(end,:,3),'o');
hold off


