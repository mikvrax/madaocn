f = dlmread('Data_Highschool2.txt');
s = size(f);
%unique time
time = unique(f(:,3));
%unique nodes
lnodes = unique(f(:,1));
rnodes = unique(f(:,2));
allnodes = [lnodes ; rnodes];
nodes = unique(allnodes);
%unique links
alllinks = f(:,1:2);
links = unique(alllinks,'rows');
%link density
p = size(links,1) / ((size(nodes,1)-1) * size(nodes,1));
%degrees,mean and variance
ED = 2*size(links,1) / size(nodes,1);
d = zeros(size(nodes,1),1);
for i=1:size(links,1)
d(links(i,1)) = d(links(i,1)) + 1;
d(links(i,2)) = d(links(i,2)) + 1;
end
vd = var(d);
md = mean(d);
%plot degree distribution
[a,b] = ksdensity(d);
figure;
plot(b,a)
A=[];
%maybe A needs A(i,i) = 1
for i=1:size(links,1)
    A(links(i,1),links(i,2)) = 1;
    A(links(i,2),links(i,1)) = 1;
end
neigh = [];
for i=1:size(nodes,1)
    k = 1;
    for j=1:size(nodes,1)
        if A(i,j) == 1
           neigh(i,k) = j;
           k = k + 1;
        end
    end
end
%find neighbours and links between them
%check if Ln correct, seems too large
Ln=zeros(1,size(nodes,1));
for i=1:size(nodes,1)
    for j=1:size(neigh(i,:),2)
        if (neigh(i,j) ~= 0)
           for k=j+1:size(neigh(i,:),2)
               if (neigh(i,k) == 0)
                   break
               end
               if (A(neigh(i,j),neigh(i,k)) == 1)
                   Ln(i) = Ln(i) + 1;
               end
           end
        else
            break
        end
    end
end
clust = [];
for i=1:size(nodes,1)
    clust(i) = Ln(i)/(d(i) * (d(i) - 1)/2);
end
clustG = 0;
for i=1:size(nodes,1)
    clustG = clustG + clust(i);
end
clustG = clustG / size(clust,2);

cd matlab_bgl-4.0.1/matlab_bgl/
SA = sparse(A);
H = all_shortest_paths(SA);
Hmax = max(max(H));
EH = mean(mean(H));


d3 = [];
for i=1:size(d,1)
    d3(i) = d(i,1) * d(i,1) * d(i,1);
end
N1 = 2*size(links,1);
N2 = d' * d;
N3 = d'*A*d;
r = (N1*N3 - N2*N2)/(N1*sum(d3) - N2*N2);

%check small-world property
ER = erdos_reyni(327,0.0546);
cER = clustering_coefficients(ER);
      
clust_ER = 0;
for i=1:size(cER,1)
clust_ER = clust_ER + cER(i);
end
clust_ER = clust_ER / size(cER,1);

HER = all_shortest_paths(ER);
EHER = mean(mean(HER));
sigma = (clustG/clust_ER)/(EH/EHER);


%eigenvalues
e = eig(A);
l1 = max(e);
D = zeros(size(nodes,1));
for i=1:size(d,1)
    D(i,i) = d(i);
end
Q = D - A;
m = eig(Q);
mn = m(size(nodes,1)-1);


%B
I=zeros(size(nodes,1),7375);
for i=1:size(nodes,1)
    t = 1;
    infe = [];
    infe(1) = i;
    k=2;
    for j=1:size(f,1)
        if f(j,3) == t
            for z=1:size(infe,2)
                if (infe(z) == f(j,1))
                    exists = ismember(f(j,2),infe);
                    if exists == 0
                        infe(k) = f(j,2);
                        k = k + 1;
                    end
                end
                if (f(j,2) == infe(z))
                    exists = ismember(f(j,1),infe);
                    if exists == 0
                        infe(k) = f(j,1);
                        k = k + 1;
                    end
                end
            end
       else
            I(i,t) = size(infe,2);
            t=t+1;
            for z=1:size(infe,2)
                if (infe(z) == f(j,1))
                    exists = ismember(f(j,2),infe);
                    if exists == 0
                        infe(k) = f(j,2);
                        k = k + 1;
                    end
                end
                if (f(j,2) == infe(z))
                    exists = ismember(f(j,1),infe);
                    if exists == 0
                        infe(k) = f(j,1);
                        k = k + 1;
                    end
                end
            end
        end
    end
    I(i,t) = size(infe,2);
end

%9
t = 1:size(time,1);
EI = mean(I);
figure;
plot(t, EI,'color','b');
e = sqrt(var(I));
hold on;
plot(t,e,'color','r');
legend("Average number of infected nodes", "Error bar")

%10
infl = 80/100 * size(nodes,1);
R= zeros(size(nodes,1),1)*size(time,1);
for i=1:size(I,1)
    for j=1:size(I,2)
        if I(i,j) >= infl
            R(i) = j;
            break;
        end
    end 
end
[Rank, IndR] = sort(R);
%11
f1 = 0.05:0.05:0.5;
[Dr, IndD] = sort(d);
[Cr, IndC] = sort(clust);
rRD = zeros(size(f1,2),1);
rRC = zeros(size(f1,2),1);
for i=1:size(f1,2)
fr = floor(f1(i) * size(Dr,1));
Drf = IndD(1:fr);
Crf = IndC(1:fr);
Rankf = IndR(1:fr);
rRD(i) = size(intersect(Rankf,Drf),1) / size(Rankf,1);
rRC(i) = size(intersect(Rankf,Crf),1) / size(Rankf,1);
end
figure;
plot(f1,rRD,'color','b')
hold on;
plot(f1,rRC,'color','r')
legend("recognition rate using degree","recognition rate using clustering coefficient")
%12
%Closeness centrality
l = zeros(size(nodes,1),1);
for i=1:size(nodes,1)
for j=1:size(nodes,1)
if (i ~= j)
l(i) = l(i) + 1/H(i,j);
end
end
end
f1 = 0.05:0.05:0.5;
[lr, Indl] = sort(l);
rRl = zeros(size(f1,2),1);
for i=1:size(f1,2)
fr = floor(f1(i) * size(lr,1));
lrf = Indl(1:fr);
Rankf = IndR(1:fr);
rRl(i) = size(intersect(Rankf,lrf),1) / size(Rankf,1);
end
figure;
plot(f1,rRl,'color','b')
legend("recognition rate using closeness centrality")
%temporal degree
td = zeros(size(nodes,1),size(time,1));
t = 1;
for i=1:size(f,1)
if (f(i,3) == t)
td(f(i,1),t) = td(f(i,1),t) + 1;
td(f(i,2),t) = td(f(i,2),t) + 1;
else
t = t+1;
td(f(i,1),t) = td(f(i,1),t) + 1;
td(f(i,2),t) = td(f(i,2),t) + 1;
end
end
etd = mean(td,2);
f1 = 0.05:0.05:0.5;
[etdr, Indetd] = sort(etd);
rRetd = zeros(size(f1,2),1);
for i=1:size(f1,2)
fr = floor(f1(i) * size(etdr,1));
etdrf = Indetd(1:fr);
Rankf = IndR(1:fr);
rRetd(i) = size(intersect(Rankf,etdrf),1) / size(Rankf,1);
end
figure;
plot(f1,rRetd,'color','b')
legend("recognition rate using mean temporal degree")

%14
g2 = f;
for i=1:size(f,1)
    ri = floor(rand() * size(f,1)) + 1;
    if (ri ~= i)
        swap = g2(i,3);
        g2(i,3) = g2(ri,3);
        g2(ri,3) = swap;
    end
end


g3 = zeros(size(f,1),3);
for i=1:size(f,1)
    ri = floor(rand() * size(f,1)) + 1;
    ri2 = floor(rand() * size(links,1)) + 1;
    g3(i,:) = [links(ri2,:) f(ri,3)];
end





g2 = sortrows(g2,3);
g3 = sortrows(g3,3);

I2=zeros(size(nodes,1),7375);
for i=1:size(nodes,1)
    t = 1;
    infe = [];
    infe(1) = i;
    k=2;
    for j=1:size(g2,1)
        if g2(j,3) == t
            for z=1:size(infe,2)
                if (infe(z) == g2(j,1))
                    exists = ismember(g2(j,2),infe);
                    if exists == 0
                        infe(k) = g2(j,2);
                        k = k + 1;
                    end
                end
                if (g2(j,2) == infe(z))
                    exists = ismember(g2(j,1),infe);
                    if exists == 0
                        infe(k) = g2(j,1);
                        k = k + 1;
                    end
                end
            end
       else
            I2(i,t) = size(infe,2);
            t=t+1;
            for z=1:size(infe,2)
                if (infe(z) == g2(j,1))
                    exists = ismember(g2(j,2),infe);
                    if exists == 0
                        infe(k) = g2(j,2);
                        k = k + 1;
                    end
                end
                if (g2(j,2) == infe(z))
                    exists = ismember(g2(j,1),infe);
                    if exists == 0
                        infe(k) = g2(j,1);
                        k = k + 1;
                    end
                end
            end
        end
    end
    I2(i,t) = size(infe,2);
end


t = 1:size(time,1);
EI2 = mean(I2);
figure;

plot(t, EI2,'color','b');
e2 = sqrt(var(I2));
hold on;
plot(t,e2,'color','r');
legend("Average number of infected nodes", "Error bar")

I3=zeros(size(nodes,1),7375);
for i=1:size(nodes,1)
    t = 1;
    infe = [];
    infe(1) = i;
    k=2;
    for j=1:size(g3,1)
        if g3(j,3) == t
            for z=1:size(infe,2)
                if (infe(z) == g3(j,1))
                    exists = ismember(g3(j,2),infe);
                    if exists == 0
                        infe(k) = g3(j,2);
                        k = k + 1;
                    end
                end
                if (g3(j,2) == infe(z))
                    exists = ismember(g3(j,1),infe);
                    if exists == 0
                        infe(k) = g3(j,1);
                        k = k + 1;
                    end
                end
            end
       else
            I3(i,t) = size(infe,2);
            t=t+1;
            for z=1:size(infe,2)
                if (infe(z) == g3(j,1))
                    exists = ismember(g3(j,2),infe);
                    if exists == 0
                        infe(k) = g3(j,2);
                        k = k + 1;
                    end
                end
                if (g3(j,2) == infe(z))
                    exists = ismember(g3(j,1),infe);
                    if exists == 0
                        infe(k) = g3(j,1);
                        k = k + 1;
                    end
                end
            end
        end
    end
    I3(i,t) = size(infe,2);
end

t = 1:size(time,1);
EI3 = mean(I3);
figure;

plot(t, EI3,'color','b');
e3 = sqrt(var(I3));
hold on;
plot(t,e3,'color','r');
legend("Average number of infected nodes", "Error bar")

