%read dataset into matrix of 4 columns
%column 1 has the first node of a link
%column 2 has the second node of a link
%column 3 has the rating of a link
%column 4 has the timestamp of a link
f = csvread("soc-sign-bitcoinotc.csv");

%create adjacency matrix
A = zeros(6005,6005);

%set the timestep of the network to one week in seconds
timestep = 604800;

%initialize matrices that stores the number of well-balanced and weak-balanced triangles 
%each node belongs to in each timestep accordingly
baln = zeros(6005,size(f,1));
unbaln = zeros(6005,size(f,1));

%initialize matrices that store the number of well-balanced and
%weak-balanced triangles in each timestep for the whole network
bal = zeros(size(f,1),1); 
unbal = zeros(size(f,1),1);

%initialize average rating of each node in well-balanced and weak-balanced triangles he is part of per week
avg_baln = zeros(6005,size(f,1));
avg_unbaln = zeros(6005,size(f,1));

%initialize average number of well and weak-balanced triangles per week
avg_bal = zeros(size(f,1),1);
avg_unbal = zeros(size(f,1),1);

stats = [];
w = 1;
nodes = [];
id = 1;
%variable that keeps track of the last processed edge between consecutive weeks
p=1;
%get the minimum and maximum timestamps
mintime = min(f(:,4));
maxtime = max(f(:,4));
%initialize matrix that holds the triangles and the number of triangles
triangles = []; 
%index where the current triangle will be stored in the matrix
l = 1;
%find the number of weeks between the minimum and the maximum timestamp
maxstep = round((maxtime-mintime)/timestep);
%end of 2013-1388527200
%beginning of 2012-1325368800
%repeat calculation of triangles for every week
for t=1:7899
    A(f(t,1),f(t,2)) = f(t,3);
    % keep track of nodes that rated each other
    nodes(1,id) = f(t,1);
    id = id + 1;
    nodes(1,id) = f(t,2);
    id = id + 1;
end

nodes = sort(unique(nodes));
id = size(nodes,2) + 1;
p=7900;

%for t=1:maxstep
for t=1:104
    %fill adjacency matrix until week t
    for i=p:size(f,1)
        %if f(i,4) <= mintime + t* timestep
        if f(i,4) <= 1325368800 + t* timestep
            A(f(i,1),f(i,2)) = f(i,3);
            % keep track of nodes that rated each other
            nodes(1,id) = f(i,1);
            id = id + 1;
            nodes(1,id) = f(i,2);
            id = id + 1;
        else
            p = i;
            break
        end
    end

    %get the number of nodes in the adjacency matrix
    s = size(A);
    nodes = sort(unique(nodes));
    id = size(nodes,2) + 1;
    %variable that keeps track of last triangle glimpsed within the previous week
    prev = l;
    %iterate over all nodes in the graph to find a third node that has
    %rated or has been rated by nodes j or k
    for u=1:id-1
        m = nodes(1,u);
        for o=1:id-1 
            j = nodes(1,o);
            %if there is a rating between nodes m and j
            if (A(m,j) ~= 0)
                for e=1:id-1
                    k = nodes(1,e);
                    %if there is a rating between nodes k and m
                    if (A(k,m) ~= 0) && (A(j,k) ~= 0)
                        %store the nodes that form the triangle, the
                        %average of their weight, the sign of the
                        %triangle and the week number               
                        [triangles(l,1),triangles(l,2),triangles(l,3)] = sorttriangle(m,j,k);
                        triangles(l,4) = (A(m,j) + A(j,k) + A(k,m))/3;
                        triangles(l,5) = sign(A(m,j)) * sign(A(j,k)) * sign(A(k,m));
                        %update index for the next triangle
                        l = l + 1;
                    end
                    %if there is a rating between nodes m and k
                    if (A(m,k) ~= 0) && (A(k,j) ~= 0)
                        %store the nodes that form the triangle, the
                        %average of their weight, the sign of the
                        %triangle and the week number
                        [triangles(l,1),triangles(l,2),triangles(l,3)] = sorttriangle(m,j,k);
                        triangles(l,4) = (A(m,j) + A(k,j) + A(m,k))/3;
                        triangles(l,5) = sign(A(m,j)) * sign(A(k,j)) * sign(A(m,k));
                        %update index for the next triangle
                        l = l + 1;
                    end
                    %if there is a rating between nodes k and m
                    if (A(k,m) ~= 0) && (A(k,j) ~= 0)
                        %store the nodes that form the triangle, the
                        %average of their weight, the sign of the
                        %triangle and the week number               
                        [triangles(l,1),triangles(l,2),triangles(l,3)] = sorttriangle(m,j,k);
                        triangles(l,4) = (A(m,j) + A(k,j) + A(k,m))/3;
                        triangles(l,5) = sign(A(m,j)) * sign(A(k,j)) * sign(A(k,m));
                        %update index for the next triangle
                        l = l + 1;
                    end
                    %if there is a rating between nodes m and k
                    if (A(m,k) ~= 0) && (A(k,j) ~= 0)
                        %store the nodes that form the triangle, the
                        %average of their weight, the sign of the
                        %triangle and the week number
                        [triangles(l,1),triangles(l,2),triangles(l,3)] = sorttriangle(m,j,k);
                        triangles(l,4) = (A(m,j) + A(k,j) + A(m,k))/3;
                        triangles(l,5) = sign(A(m,j)) * sign(A(k,j)) * sign(A(m,k));
                        %update index for the next triangle
                        l = l + 1;
                    end
                end
            end
            %if there is a rating between nodes j and m
            if (A(j,m) ~= 0)
                for e=1:id-1
                    k = nodes(1,e);
                    %if there is a rating between nodes k and m
                    if (A(k,m) ~= 0) && (A(j,k) ~= 0)  
                        %store the nodes that form the triangle, the
                        %average of their weight, the sign of the
                        %triangle and the week number
                        [triangles(l,1),triangles(l,2),triangles(l,3)] = sorttriangle(m,j,k);
                        triangles(l,4) = (A(j,m) + A(j,k) + A(k,m))/3;
                        triangles(l,5) = sign(A(j,m)) * sign(A(j,k)) * sign(A(k,m));
                        %update index for the next triangle
                        l = l + 1;
                    end
                    %if there is a rating between nodes m and k
                    if (A(m,k) ~= 0) && (A(j,k) ~= 0) 
                        %store the nodes that form the triangle, the
                        %average of their weight, the sign of the
                        %triangle and the week number
                        [triangles(l,1),triangles(l,2),triangles(l,3)] = sorttriangle(m,j,k);
                        triangles(l,4) = (A(j,m) + A(j,k) + A(m,k))/3;
                        triangles(l,5) = sign(A(j,m)) * sign(A(j,k)) * sign(A(m,k));
                        %update index for the next triangle
                        l = l + 1;
                    end
                    %if there is a rating between nodes k and m
                    if (A(k,m) ~= 0) && (A(k,j) ~= 0)  
                        %store the nodes that form the triangle, the
                        %average of their weight, the sign of the
                        %triangle and the week number
                        [triangles(l,1),triangles(l,2),triangles(l,3)] = sorttriangle(m,j,k);
                        triangles(l,4) = (A(j,m) + A(k,j) + A(k,m))/3;
                        triangles(l,5) = sign(A(j,m)) * sign(A(k,j)) * sign(A(k,m));
                        %update index for the next triangle
                        l = l + 1;
                    end
                    %if there is a rating between nodes m and k
                    if (A(m,k) ~= 0) && (A(k,j) ~= 0) 
                        %store the nodes that form the triangle, the
                        %average of their weight, the sign of the
                        %triangle and the week number
                        [triangles(l,1),triangles(l,2),triangles(l,3)] = sorttriangle(m,j,k);
                        triangles(l,4) = (A(j,m) + A(k,j) + A(m,k))/3;
                        triangles(l,5) = sign(A(j,m)) * sign(A(k,j)) * sign(A(m,k));
                        %update index for the next triangle
                        l = l + 1;
                    end
                end
            end
        end
    end
    if (size(triangles,1) > 0)
        utriangles = unique(triangles,'rows');
        triangles = [];
        triangles(1,:) = utriangles(1,:);
        q = 2;
        for i=2:size(utriangles,1)
            if ((utriangles(i,1) == utriangles(i-1,1)) && (utriangles(i,2) == utriangles(i-1,2)) && (utriangles(i,3) == utriangles(i-1,3)))
                if (utriangles(i,5) == utriangles(i-1,5))
                    triangles(q-1,4) = (triangles(q-1,4) + utriangles(i,4))/2;
                else
                    triangles(q,:) = utriangles(i,:);
                    q = q + 1;
                end
            else
                triangles(q,:) = utriangles(i,:);
                q = q + 1;
            end
        end
        l = q;
        % get unique triangles(combine triangles that share the same nodes, rating and type)

        % calculate over time the number of balanced and unbalanced triangles in the network,
        % the average rating of the balanced and the unbalanced triangles
        % the same calculations are done for each node over time
        for i=1:size(triangles,1)
            ts = t;
            m = triangles(i,1);
            j = triangles(i,2);
            k = triangles(i,3);
            if triangles(i,5) == 1  
                 avg_baln(m,ts) = avg_baln(m,ts) + A(m,j);
                 avg_baln(j,ts) = avg_baln(j,ts) + A(j,k);
                 avg_baln(k,ts) = avg_baln(k,ts) + A(k,m);
                 avg_baln(k,ts) = avg_baln(k,ts) + A(m,k);
                 baln(m,ts) = baln(m,ts) + 1;
                 baln(j,ts) = baln(j,ts) + 1;
                 baln(k,ts) = baln(k,ts) + 1;
            else
                 avg_unbal(m,ts) = avg_unbaln(m,ts) + A(m,j);
                 avg_unbaln(j,ts) = avg_unbaln(j,ts) + A(j,k);
                 avg_unbaln(k,ts) = avg_unbaln(k,ts) + A(k,m);
                 unbaln(m,ts) = unbaln(m,ts) + 1;
                 unbaln(j,ts) = unbaln(j,ts) + 1;
                 unbaln(k,ts) = unbaln(k,ts) + 1;
            end
        end

        
        for i=1:size(nodes,2)
            stats(w,1) = nodes(1,i);
            stats(w,2) = baln(nodes(1,i),t);
            stats(w,3) = unbaln(nodes(1,i),t);
            %calculate average rating of nodes among the triangles they are part of over time
            if (baln(nodes(1,i),t) ~= 0)
                stats(w,4) = avg_baln(nodes(1,i),t)/baln(nodes(1,i),t);
            else
                stats(w,4) = 0;
            end
            if (unbaln(nodes(1,i),t) ~= 0)
                stats(w,5) = avg_unbaln(nodes(1,i),t)/unbaln(nodes(1,i),t);
            else
                stats(w,5) = 0;
            end
            if (unbaln(nodes(1,i),t) ~= 0) && (baln(nodes(1,i),t) ~= 0)
                stats(w,6) = baln(nodes(1,i),t)/(baln(nodes(1,i),t) + unbaln(nodes(1,i),t));
            else
                stats(w,6) = 0;
            end
            w = w + 1;
        end
    end
end
csvwrite('triangles_stats.csv',stats);
save('triangles3.mat');
