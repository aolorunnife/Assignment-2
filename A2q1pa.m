L = 600E-9;
W = 600E-9;
vo = 1;

nx = 30;
ny = 30;

V = zeros(nx,ny);
F = zeros(nx*ny,1);
G = sparse(nx*ny,nx*ny);

for j = 1:ny
    for i = 1:nx

        n = j+(i-1)*ny;
        if i == 1
            F(n) = vo;

        elseif i == nx
            F(n) = 0;

        end

        if i == nx || i == 1
            G(n,n) = 1;


        elseif j == ny
            nym = j-1+(i-1)*ny;
            nxm = j+(i-1-1)*ny;
            nxp = j+(i-1+1)*ny;
            G(n,n) = -3;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;

        elseif j == 1
            nxp = j+(i-1+1)*ny;
            nxm = j+(i-1-1)*ny;
            nyp = j+1+(i-1)*ny;
            G(n,n) = -3;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nyp) = 1;

        else
            nyp = j+1+(i-1)*ny;
            nxm = j+(i-1-1)*ny;
            nxp = j+(i-1+1)*ny;
            nym = j-1+(i-1)*ny;
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nyp) = 1;
            G(n,nym) = 1;

        end
    end

end


for j = 1:ny
    for i=1:nx
        P = G\F;
        n = j+(i-1)*ny;
        V(i,j) = P(n);
    end
end

axis([0 L 0 W])
surf(V)
xlabel('x');
ylabel('y');
zlabel('z');
title(sprintf('Finite Differences (1D) Solution in (3D)'));


[Ex, Ey] = gradient(V);
figure, quiver(-Ey', -Ex', 10)
title(sprintf('Gradient Field Solution'));
xlabel('x');
ylabel('y');

