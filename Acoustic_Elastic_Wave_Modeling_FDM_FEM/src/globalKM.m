function [ K, M ] = globalKM( nz, Ke, Me )

% Global matrices have a dimension with (nz^2 x nz^2).
nz = floor(nz);
N = floor(nz * nz);
K = zeros(N, N);
M = zeros(N, N);

% Assemble global stiffness and mass matrices.
for i = 1: N
    % 1. Lower left corner.
   if i == 1
       K(i, i) = Ke(1, 1); K(i, i+1) = Ke(1, 4); K(i, i+nz) = Ke(1, 2); K(i, i+nz+1) = Ke(1, 3);
       M(i, i) = Me(1, 1); M(i, i+1) = Me(1, 4); M(i, i+nz) = Me(1, 2); M(i, i+nz+1) = Me(1, 3);
   % 2. Lower right corner.
   elseif i == nz
       K(i, i-1) = Ke(4, 1); K(i, i) = Ke(4, 4); K(i, i+nz-1) = Ke(4, 2); K(i, i+nz) = Ke(4, 3);
       M(i, i-1) = Me(4, 1); M(i, i) = Me(4, 4); M(i, i+nz-1) = Me(4, 2); M(i, i+nz) = Me(4, 3);
   % 3. Upper left corner.
   elseif i == ((nz-1)*nz+1)
       K(i, i-nz) = Ke(2, 1); K(i, i-nz+1) = Ke(2, 4); K(i, i) = Ke(2, 2); K(i, i+1) = Ke(2, 3);
       M(i, i-nz) = Me(2, 1); M(i, i-nz+1) = Me(2, 4); M(i, i) = Me(2, 2); M(i, i+1) = Me(2, 3);
   % 4. Upper right corner.
   elseif i == nz*nz
       K(i, i-nz-1) = Ke(3, 1); K(i, i-nz) = Ke(3, 4); K(i, i-1) = Ke(3, 2); K(i, i) = Ke(3, 3);
       M(i, i-nz-1) = Me(3, 1); M(i, i-nz) = Me(3, 4); M(i, i-1) = Me(3, 2); M(i, i) = Me(3, 3);
   % 5. Left column.
   elseif i > nz && i < (nz-1)*nz && mod((i-1),nz)==0
       K(i, i-nz) = Ke(2, 1); K(i, i-nz+1) = Ke(2, 4); K(i, i) = Ke(1, 1) + Ke(2, 2);
       K(i, i+1) = Ke(2, 3) + Ke(1, 4); K(i, i+nz) = Ke(1, 2); K(i, i+nz+1) = Ke(1, 3);
       
       M(i, i-nz) = Me(2, 1); M(i, i-nz+1) = Me(2, 4); M(i, i) = Me(1, 1) + Me(2, 2);
       M(i, i+1) = Me(2, 3) + Me(1, 4); M(i, i+nz) = Me(1, 2); M(i, i+nz+1) = Me(1, 3);
   % 6. Right column.
   elseif i > nz && i < N && mod(i, nz)==0
       K(i, i-nz-1) = Ke(3, 1); K(i, i-nz) = Ke(3, 4); K(i, i-1) = Ke(4, 1) + Ke(3, 2); 
       K(i, i) = Ke(3, 3) + Ke(4, 4); K(i, i+nz-1) = Ke(4, 2); K(i, i+nz) = Ke(4, 3);
       
       M(i, i-nz-1) = Me(3, 1); M(i, i-nz) = Me(3, 4); M(i, i-1) = Me(4, 1) + Me(3, 2); 
       M(i, i) = Me(3, 3) + Me(4, 4); M(i, i+nz-1) = Me(4, 2); M(i, i+nz) = Me(4, 3);
   % 7. Lower row.
   elseif i > 1 && i < nz
       K(i, i-1) = Ke(4, 1); K(i, i) = Ke(1, 1) + Ke(4, 4); K(i, i+1) = Ke(1, 4); 
       K(i, i+nz-1) = Ke(4, 2); K(i, i+nz) = Ke(1, 2) + Ke(4, 3); K(i, i+nz+1) = Ke(1, 3);
       
       M(i, i-1) = Me(4, 1); M(i, i) = Me(1, 1) + Me(4, 4); M(i, i+1) = Me(1, 4); 
       M(i, i+nz-1) = Me(4, 2); M(i, i+nz) = Me(1, 2) + Me(4, 3); M(i, i+nz+1) = Me(1, 3);
       
   % 8. Upper row.
   elseif i > ((nz-1)*nz+1) && i < N
       K(i, i-nz-1) = Ke(3, 1); K(i, i-nz) = Ke(2, 1) + Ke(3, 4); K(i, i-nz+1) = Ke(2, 4); 
       K(i, i-1) = Ke(3, 2); K(i, i) = Ke(2, 2) + Ke(3, 3); K(i, i+1) = Ke(2, 3);
       
       M(i, i-nz-1) = Me(3, 1); M(i, i-nz) = Me(2, 1) + Me(3, 4); M(i, i-nz+1) = Me(2, 4); 
       M(i, i-1) = Me(3, 2); M(i, i) = Me(2, 2) + Me(3, 3); M(i, i+1) = Me(2, 3);
   % 9. Internal part.
   else
       K(i, i-nz-1) = Ke(3, 1); K(i, i-nz) = Ke(2, 1) + Ke(3, 4); K(i, i-nz+1) = Ke(2, 4);
       K(i, i-1) = Ke(4, 1) + Ke(3, 2); K(i, i) = Ke(1, 1) + Ke(2, 2) + Ke(3, 3) + Ke(4, 4); 
       K(i, i+1) = Ke(1, 4) + Ke(2, 3); K(i, i+nz-1) = Ke(4, 2); 
       K(i, i+nz) = Ke(1, 2) + Ke(4, 3); K(i, i+nz+1) = Ke(1, 3);
       
       M(i, i-nz-1) = Me(3, 1); M(i, i-nz) = Me(2, 1) + Me(3, 4); M(i, i-nz+1) = Me(2, 4);
       M(i, i-1) = Me(4, 1) + Me(3, 2); M(i, i) = Me(1, 1) + Me(2, 2) + Me(3, 3) + Me(4, 4); 
       M(i, i+1) = Me(1, 4) + Me(2, 3); M(i, i+nz-1) = Me(4, 2); 
       M(i, i+nz) = Me(1, 2) + Me(4, 3); M(i, i+nz+1) = Me(1, 3);
   end
end
end