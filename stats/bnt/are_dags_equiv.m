function equiv = are_dags_equiv( G1, G2)
% Check if two graphs are equivalent. Use Pearl & Verma criteria,
% i.e. check for same skeleton and V-structure


% Check for skeleton
equiv = isequal(G1+G1', G2+G2');

if(equiv==0)
    return;
end


n = size(G1,1);
% Check for V-structure
for i=1:n
    for j=setdiff(1:n,i)
        for k=setdiff(1:n, [i,j])

            if( (G1(i,k)*G1(j,k)* (1-G1(i,j))*(1-G1(j,i))) &&   (1-(G2(i,k)*G2(j,k)* (1-G2(i,j))*(1-G2(j,i))))   ) % Check V-structure in G1
                equiv=0;  return;
            end

            if( (G2(i,k)*G2(j,k)* (1-G2(i,j))*(1-G2(j,i))) &&   (1-(G1(i,k)*G1(j,k)* (1-G1(i,j))*(1-G1(j,i))))   ) % Check V-structure in G2
                equiv=0;  return;
            end
        end
    end
end


