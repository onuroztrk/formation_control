V = [-1 -1 -1;1 -1 -1;1 1 -1;-1 1 -1;-1 -1 1;1 -1 1;1 1 1;-1 1 1;0 0 0];
d14 = 2;
d15 = 2;
d19 = 1.7321;
d23 = 2;
d26 = 2;
d29 = 1.7321;
d34 = 2;
d37 = 2;
d39 = 1.7321;
d45 = 2*sqrt(2);
d48 = 2;
d49 = 1.7321;
d56 = 2;
d58 = 2;
d59 = 1.7321;
d67 = 2;
d68 = 2*sqrt(2);
d69 = 1.7321;
d78 = 2;
d79 = 1.7321;
d89 = 1.7321;
dij_ = zeros(1,21);
ddij_ = zeros(1,21);
kv = 1;
ka = 1;
T = 0.01;
q0 = zeros(9,3);
v0 = zeros(9,3);
qout = zeros(3000,9,3);
qij_out = zeros(3000,1,21);
e_out = zeros(3000,21);
for i = 1:9
    q0(i,:) = V(i,:)+ 0.5*(rand(1,3)-0.5*ones(1,3));
    v0(i,:) = 2*(rand(1,3)-0.5*ones(1,3));
end

q = q0;
q_ = q0;
vi = v0;
qin_ = [(q0(1,:)-q0(9,:));(q0(2,:)-q0(9,:));(q0(3,:)-q0(9,:));(q0(4,:)-q0(9,:));(q0(5,:)-q0(9,:)); ...
           (q0(6,:)-q0(9,:));(q0(7,:)-q0(9,:));(q0(8,:)-q0(9,:));0 0 0];
       
Rq_ = [(q(1,:)-q(4,:)) zeros(1,6) (q(4,:)-q(1,:)) zeros(1,15);
          (q(1,:)-q(5,:)) zeros(1,9) (q(5,:)-q(1,:)) zeros(1,12);
          (q(1,:)-q(9,:)) zeros(1,21) (q(9,:)-q(1,:));
          zeros(1,3) (q(2,:)-q(3,:)) (q(3,:)-q(2,:)) zeros(1,18);
          zeros(1,3) (q(2,:)-q(6,:)) zeros(1,9) (q(6,:)-q(2,:)) zeros(1,9);
          zeros(1,3) (q(2,:)-q(9,:)) zeros(1,18) (q(9,:)-q(2,:));
          zeros(1,6) (q(3,:)-q(4,:)) (q(4,:)-q(3,:)) zeros(1,15);
          zeros(1,6) (q(3,:)-q(7,:)) zeros(1,9) (q(7,:)-q(3,:)) zeros(1,6);
          zeros(1,6) (q(3,:)-q(9,:)) zeros(1,15) (q(9,:)-q(3,:));
          zeros(1,9) (q(4,:)-q(5,:)) (q(5,:)-q(4,:)) zeros(1,12);
          zeros(1,9) (q(4,:)-q(8,:)) zeros(1,9) (q(8,:)-q(4,:)) zeros(1,3);
          zeros(1,9) (q(4,:)-q(9,:)) zeros(1,12) (q(9,:)-q(4,:));
          zeros(1,12) (q(5,:)-q(6,:)) (q(6,:)-q(5,:)) zeros(1,9);
          zeros(1,12) (q(5,:)-q(8,:)) zeros(1,6) (q(8,:)-q(5,:)) zeros(1,3);
          zeros(1,12) (q(5,:)-q(9,:)) zeros(1,9) (q(9,:)-q(5,:));
          zeros(1,15) (q(6,:)-q(7,:)) (q(7,:)-q(6,:)) zeros(1,6);
          zeros(1,15) (q(6,:)-q(8,:)) zeros(1,3) (q(8,:)-q(6,:)) zeros(1,3);
          zeros(1,15) (q(6,:)-q(9,:)) zeros(1,6) (q(9,:)-q(6,:));
          zeros(1,18) (q(7,:)-q(8,:)) (q(8,:)-q(7,:)) zeros(1,3);
          zeros(1,18) (q(7,:)-q(9,:)) zeros(1,3) (q(9,:)-q(7,:));
          zeros(1,21) (q(8,:)-q(9,:)) (q(9,:)-q(8,:))];
Rq_plus_ = Rq_'/(Rq_*Rq_');
for t = 1:3000
    
    v = [1 cos(t*T) 0];
    vdot = [0 -sin(t*T) 0];
    w = [1 1 1];
    pt = 1+0.5*sin(0.4*t*T);
    dij = pt*[d14 d15 d19 d23 d26 d29 d34 d37 d39 d45 d48 d49 d56 d58 d59 d67 d68 d69 d78 d79 d89];

    if t == 1
        dij_ = dij;
    end
    ddij = (dij-dij_)/T;
    dij_ = dij;
    if t == 1
        ddij_ = ddij;
    end
    dddij = (ddij-ddij_)/T;
    ddij_ = ddij;
    dv = dij.*ddij;
    ddv = ddij.*ddij + dij.*dddij;

    Rq = [(q(1,:)-q(4,:)) zeros(1,6) (q(4,:)-q(1,:)) zeros(1,15);
          (q(1,:)-q(5,:)) zeros(1,9) (q(5,:)-q(1,:)) zeros(1,12);
          (q(1,:)-q(9,:)) zeros(1,21) (q(9,:)-q(1,:));
          zeros(1,3) (q(2,:)-q(3,:)) (q(3,:)-q(2,:)) zeros(1,18);
          zeros(1,3) (q(2,:)-q(6,:)) zeros(1,9) (q(6,:)-q(2,:)) zeros(1,9);
          zeros(1,3) (q(2,:)-q(9,:)) zeros(1,18) (q(9,:)-q(2,:));
          zeros(1,6) (q(3,:)-q(4,:)) (q(4,:)-q(3,:)) zeros(1,15);
          zeros(1,6) (q(3,:)-q(7,:)) zeros(1,9) (q(7,:)-q(3,:)) zeros(1,6);
          zeros(1,6) (q(3,:)-q(9,:)) zeros(1,15) (q(9,:)-q(3,:));
          zeros(1,9) (q(4,:)-q(5,:)) (q(5,:)-q(4,:)) zeros(1,12);
          zeros(1,9) (q(4,:)-q(8,:)) zeros(1,9) (q(8,:)-q(4,:)) zeros(1,3);
          zeros(1,9) (q(4,:)-q(9,:)) zeros(1,12) (q(9,:)-q(4,:));
          zeros(1,12) (q(5,:)-q(6,:)) (q(6,:)-q(5,:)) zeros(1,9);
          zeros(1,12) (q(5,:)-q(8,:)) zeros(1,6) (q(8,:)-q(5,:)) zeros(1,3);
          zeros(1,12) (q(5,:)-q(9,:)) zeros(1,9) (q(9,:)-q(5,:));
          zeros(1,15) (q(6,:)-q(7,:)) (q(7,:)-q(6,:)) zeros(1,6);
          zeros(1,15) (q(6,:)-q(8,:)) zeros(1,3) (q(8,:)-q(6,:)) zeros(1,3);
          zeros(1,15) (q(6,:)-q(9,:)) zeros(1,6) (q(9,:)-q(6,:));
          zeros(1,18) (q(7,:)-q(8,:)) (q(8,:)-q(7,:)) zeros(1,3);
          zeros(1,18) (q(7,:)-q(9,:)) zeros(1,3) (q(9,:)-q(7,:));
          zeros(1,21) (q(8,:)-q(9,:)) (q(9,:)-q(8,:))];
    
    Rq_plus = Rq'/(Rq*Rq');
    dRq_plus = (Rq_plus - Rq_plus_)/T;
    Rq_plus_ = Rq_plus;

    qij = [norm((q(1,:)-q(4,:))) norm((q(1,:)-q(5,:))) norm((q(1,:)-q(9,:))) ...
           norm((q(2,:)-q(3,:))) norm((q(2,:)-q(6,:))) norm((q(2,:)-q(9,:))) ...
           norm((q(3,:)-q(4,:))) norm((q(3,:)-q(7,:))) norm((q(3,:)-q(9,:))) ...
           norm((q(4,:)-q(5,:))) norm((q(4,:)-q(8,:))) norm((q(4,:)-q(9,:))) ...
           norm((q(5,:)-q(6,:))) norm((q(5,:)-q(8,:))) norm((q(5,:)-q(9,:))) ...
           norm((q(6,:)-q(7,:))) norm((q(6,:)-q(8,:))) norm((q(6,:)-q(9,:))) ...
           norm((q(7,:)-q(8,:))) norm((q(7,:)-q(9,:))) norm((q(8,:)-q(9,:)))];
       
    qin = [(q(1,:)-q(9,:));(q(2,:)-q(9,:));(q(3,:)-q(9,:));(q(4,:)-q(9,:));(q(5,:)-q(9,:)); ...
           (q(6,:)-q(9,:));(q(7,:)-q(9,:));(q(8,:)-q(9,:));0 0 0];
    
    vin = (qin-qin_)/T;
    qin_ = qin;
       
    e = qij - dij;
    e_out(t,:) = e;
    
    z = e.*(e+2*dij);
    vf = -kv*rproduct(Rq',z) + v.*ones(9,3) + cross(w.*ones(9,3),qin);
    dz = 2*(rproduct2(Rq,vi)'-dv);
    dvf = rproduct(dRq_plus,-kv*z+dv)+rproduct(Rq_plus,-kv*dz+ddv) + vdot.*ones(9,3) + cross(w.*ones(9,3),vin);
    u = -ka*(vi-vf)+dvf-rproduct(Rq',z);

    vi = vi+u*T;
    q = q+vi*T;
    qout(t,:,:) = q;
end

% figure(1)
% plot3(qout(1,:,1),qout(1,:,2),qout(1,:,3),'.');
% hold on
% grid on
% for t = 2:1000
%     plot3(qout(t,:,1),qout(t,:,2),qout(t,:,3),'.');
% end
% plot3(qout(end,:,1),qout(end,:,2),qout(end,:,3),'o');
% hold off
g1 = graph([1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 8],[4 5 9 3 6 9 4 7 9 5 8 9 6 8 9 7 8 9 8 9 9],squeeze(qij_out(1,:,:)));
p = plot(g1,'XData',qout(1,:,1)','YData',qout(1,:,2)','ZData',qout(1,:,3)');
p.Marker = 'o';
p.NodeColor = 'r';
hold on
grid on

for i = 300:300:3000
    g1 = graph([1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 8],[4 5 9 3 6 9 4 7 9 5 8 9 6 8 9 7 8 9 8 9 9],squeeze(qij_out(i,:,:)));
    plot(g1,'XData',qout(i,:,1)','YData',qout(i,:,2)','ZData',qout(i,:,3)');
end
p.Marker = '*';
p.NodeColor = 'b';
hold off
figure(2)
plot(0:0.01:4.99,e_out(1:500,:));

function m = rproduct(R,z)
    m = zeros(9,3);
    for i = 1:9
        sum = 0;
        for j = 1:21
            sum = sum + R(3*(i-1)+1:3*(i-1)+3,j)*z(j);
        end
    m(i,:) = sum;
    end
end

function m = rproduct2(R,z)
    m = zeros(21,1);
    for i = 1:21
        sum = 0;
        for j = 1:9
            sum = sum + R(i,3*(j-1)+1:3*(j-1)+3)*z(j,:)';
        end
    m(i,:) = sum;
    end
end