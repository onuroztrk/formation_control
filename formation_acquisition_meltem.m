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
    q0(i,:) = V(i,:)+ 0.5*(rand(1,3)-0.5*ones(1,3));
    v0(i,:) = 2*(rand(1,3)-0.5*ones(1,3));
end
q = q0;
q_ = q0;
vi = v0;
Rq_ = [(q(1,:)-q(2,:)) (q(2,:)-q(1,:)) zeros(1,18);
          (q(1,:)-q(3,:)) zeros(1,3) (q(3,:)-q(1,:)) zeros(1,15);
          (q(1,:)-q(4,:)) zeros(1,6) (q(4,:)-q(1,:)) zeros(1,12);
          (q(1,:)-q(5,:)) zeros(1,9) (q(5,:)-q(1,:)) zeros(1,9);
          (q(1,:)-q(8,:)) zeros(1,18) (q(8,:)-q(1,:));
          zeros(1,3) (q(2,:)-q(3,:)) (q(3,:)-q(2,:)) zeros(1,15);
          zeros(1,3) (q(2,:)-q(5,:)) zeros(1,6) (q(5,:)-q(2,:)) zeros(1,9);
          zeros(1,3) (q(2,:)-q(6,:)) zeros(1,9) (q(6,:)-q(2,:)) zeros(1,6);
          zeros(1,6) (q(3,:)-q(4,:)) (q(4,:)-q(3,:)) zeros(1,12);
          zeros(1,6) (q(3,:)-q(6,:)) zeros(1,6) (q(6,:)-q(3,:)) zeros(1,6);
          zeros(1,6) (q(3,:)-q(7,:)) zeros(1,9) (q(7,:)-q(3,:)) zeros(1,3);
          zeros(1,9) (q(4,:)-q(7,:)) zeros(1,6) (q(7,:)-q(4,:)) zeros(1,3);
          zeros(1,9) (q(4,:)-q(8,:)) zeros(1,9) (q(8,:)-q(4,:));
          zeros(1,12) (q(5,:)-q(6,:)) (q(6,:)-q(5,:)) zeros(1,6);
          zeros(1,12) (q(5,:)-q(8,:)) zeros(1,6) (q(8,:)-q(5,:));
          zeros(1,15) (q(6,:)-q(7,:)) (q(7,:)-q(6,:)) zeros(1,3);
          zeros(1,15) (q(6,:)-q(8,:)) zeros(1,3) (q(8,:)-q(6,:));
          zeros(1,18) (q(7,:)-q(8,:)) (q(8,:)-q(7,:))];
      
for t = 1:1000

    Rq = [(q(1,:)-q(2,:)) (q(2,:)-q(1,:)) zeros(1,18);
          (q(1,:)-q(3,:)) zeros(1,3) (q(3,:)-q(1,:)) zeros(1,15);
          (q(1,:)-q(4,:)) zeros(1,6) (q(4,:)-q(1,:)) zeros(1,12);
          (q(1,:)-q(5,:)) zeros(1,9) (q(5,:)-q(1,:)) zeros(1,9);
          (q(1,:)-q(8,:)) zeros(1,18) (q(8,:)-q(1,:));
          zeros(1,3) (q(2,:)-q(3,:)) (q(3,:)-q(2,:)) zeros(1,15);
          zeros(1,3) (q(2,:)-q(5,:)) zeros(1,6) (q(5,:)-q(2,:)) zeros(1,9);
          zeros(1,3) (q(2,:)-q(6,:)) zeros(1,9) (q(6,:)-q(2,:)) zeros(1,6);
          zeros(1,6) (q(3,:)-q(4,:)) (q(4,:)-q(3,:)) zeros(1,12);
          zeros(1,6) (q(3,:)-q(6,:)) zeros(1,6) (q(6,:)-q(3,:)) zeros(1,6);
          zeros(1,6) (q(3,:)-q(7,:)) zeros(1,9) (q(7,:)-q(3,:)) zeros(1,3);
          zeros(1,9) (q(4,:)-q(7,:)) zeros(1,6) (q(7,:)-q(4,:)) zeros(1,3);
          zeros(1,9) (q(4,:)-q(8,:)) zeros(1,9) (q(8,:)-q(4,:));
          zeros(1,12) (q(5,:)-q(6,:)) (q(6,:)-q(5,:)) zeros(1,6);
          zeros(1,12) (q(5,:)-q(8,:)) zeros(1,6) (q(8,:)-q(5,:));
          zeros(1,15) (q(6,:)-q(7,:)) (q(7,:)-q(6,:)) zeros(1,3);
          zeros(1,15) (q(6,:)-q(8,:)) zeros(1,3) (q(8,:)-q(6,:));
          zeros(1,18) (q(7,:)-q(8,:)) (q(8,:)-q(7,:))];
    
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
    vf = -kv*rproduct(Rq',z);
    dz = 2*rproduct2(Rq,vi)';
    dvf = -kv*(rproduct(vij,z)+rproduct(Rq',dz));
    u = -ka*(vi-vf)+dvf-rproduct(Rq',z);
%     for i = 1:8
%         sum = 0;
%         for j = 1:18
%              sum = sum + (ka*kv + 1)*Rqt(3*(i-1)+1:3*(i-1)+3,j)*z(j)+ (kv*(z(j) + 2*Rqt(3*(i-1)+1:3*(i-1)+3,j)*Rqt(3*(i-1)+1:3*(i-1)+3,j)')*vij(3*(i-1)+1:3*(i-1)+3,j));
%         end
%         u(i,:) = -ka*vi(i,:) - sum';
%     end

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

function m = rproduct(R,z)
    m = zeros(8,3);
    for i = 1:8
        sum = 0;
        for j = 1:18
            sum = sum + R(3*(i-1)+1:3*(i-1)+3,j)*z(j);
        end
    m(i,:) = sum;
    end
end

function m = rproduct2(R,z)
    m = zeros(18,1);
    for i = 1:18
        sum = 0;
        for j = 1:8
            sum = sum + R(i,3*(j-1)+1:3*(j-1)+3)*z(j,:)';
        end
    m(i,:) = sum;
    end
end

