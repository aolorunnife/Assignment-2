L = 600E-9; %x direction
W = 400E-9; %y direction
vo = 1;

for f = 1:12

    nx = 30*f;
    ny = 20*f;
    COut = 1; %conductivity outside boxs
    CIn = 10^-2; %conductivity inside boxs
    box = 1; %making box

    V = zeros(nx,ny);
    F = zeros(nx*ny,1);
    G = sparse(nx*ny,nx*ny);  %matrix
    CM = zeros(nx,ny);

    if box == 1
        b1 = rectangle('Position',[2.0E-7,0,1.8E-7,1.8E-7]); %box1 position (bottom)
        b2 = rectangle('Position',[2.0E-7,2.2E-7,1.8E-7,1.8E-7]); %box2 position (top)
    end


    for j = 1:ny
        for i = 1:nx
            newW = j*(W/ny);
            newL = i*(L/nx);

            if box == 1 && (newW > 2.2E-7 || newW < 1.8E-7) && box == 1 && newL > 2.0E-7 && newL < 3.95E-7
                CM(i,j) = CIn;
            else
                CM(i,j) = COut;
            end
        end
    end


    eq = 1;
    for j = 1:ny
        for i = 1:nx
            n = j + (i-1)*ny;

            if i == 1
                F(n) = vo;
                G(n,n) = CM(i,j);

            elseif i == nx
                if eq == 1
                    F(n) = 0;
                    G(n,n) = CM(i,j);
                else
                    F(n) = vo;
                    G(n,n) = CM(i,j);
                end


            elseif j == 1
                if eq == 1
                    F(n) = 0;
                    nxm = j + ((i-1)-1)*ny; %(i-1,j)
                    nyp = (j+1) + (i-1)*ny; %(i,j+1)
                    nxp = j + ((i+1)-1)*ny; %(i+1,j)
                    G(n,n) = -(CM(i-1,j) + CM(i+1,j) + CM(i,j+1));
                    G(n,nxm) = CM(i-1,j);
                    G(n,nxp) = CM(i+1,j);
                    G(n,nyp) = CM(i,j+1);
                else
                    F(n) = 0;
                    G(n,n) = CM(i,j);
                end


            elseif j == ny
                if eq == 1
                    F(n) = 0;
                    nxm = j + ((i-1)-1)*ny; %(i-1,j)
                    nxp = j + ((i+1)-1)*ny; %(i+1,j)
                    nym = (j-1) + (i-1)*ny; %(i,j-1)
                    G(n,n) = -(CM(i-1,j) + CM(i+1,j) + CM(i,j-1));
                    G(n,nxm) = CM(i-1,j);
                    G(n,nxp) = CM(i+1,j);
                    G(n,nym) = CM(i,j-1);
                else
                    F(n) = 0;
                    G(n,n) = CM(i,j);
                end

            else
                nxm = j + ((i-1)-1)*ny;
                nxp = j + ((i+1)-1)*ny;
                nym = (j-1) + (i-1)*ny;
                nyp = (j+1) + (i-1)*ny;

                G(n,n) = -(CM(i-1,j) + CM(i+1,j) + CM(i,j-1) + CM(i,j+1));
                G(n,nxm) = CM(i-1,j);
                G(n,nxp) = CM(i+1,j);
                G(n,nym) = CM(i,j-1);
                G(n,nyp) = CM(i,j+1);

            end
        end
    end

    P = G\F;

    for j = 1:ny
        for i = 1:nx
            n = j + (i-1)*ny;
            V(i,j) = P(n);
        end
    end


    [Ey, Ex] = gradient(V);

    Jx = -CM.*Ex;
    Jy = -CM.*Ey;
    depth = 1;

    %Current
    A=depth.*W;

    Ix(f)=mean(Jx(1,:)).*A;

end

%plotting current vs meshize
plot((1:12)*30, Ix)  
title('Current vs Mesh Size')
xlabel('Mesh Size (nx)')
ylabel('Current (A)')


