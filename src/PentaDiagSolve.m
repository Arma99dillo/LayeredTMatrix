function x = PentaDiagSolve(A, b, m)
%solve the linear system using the band gaussian elimination with five
%bands of blocks of dimension m

nB = size(A,1)/m; % number of blocks
A0 = cell(nB,1);  % diagonal
B  = cell(nB-1,1); % first upper diagonal
C  = cell(nB-2,1); % second upper diagonal
D  = cell(nB-1,1); % first lower diagonal
E  = cell(nB-2,1); % second lower diagonal

% Block extraction
for i = 1:nB
    ii = (i-1)*m + (1:m);
    A0{i} = full(A(ii,ii));
    if i < nB
        jj = i*m + (1:m);
        B{i} = full(A(ii,jj));
        D{i} = full(A(jj,ii));
    end
    if i < nB-1
        kk = (i+1)*m + (1:m);
        C{i} = full(A(ii,kk));
        E{i} = full(A(kk,ii));
    end
end

% RHS as blocks
bb = cell(nB,1);
for i = 1:nB
    bb{i} = b((i-1)*m + (1:m));
end

% ===============================
% FORWARD ELIMINATION
% ===============================
for k = 1:nB-1
    % --- Step 1: Eliminate D{k} (affects row k+1) ---
    % If k=1 and A0{1}=I, Mk1 is just D{1}. Otherwise, use / operator.
    if k == 1
        Mk1 = D{k};
    else
        Mk1 = D{k} / A0{k};
    end

    A0{k+1} = A0{k+1} - Mk1 * B{k};
    if k < nB-1
        B{k+1} = B{k+1} - Mk1 * C{k};
    end
    bb{k+1} = bb{k+1} - Mk1 * bb{k};

    % --- Step 2: Eliminate E{k} (affects row k+2) ---
    if k < nB-1
        if k == 1
            Mk2 = E{k};
        else
            Mk2 = E{k} / A0{k};
        end

        % This is the critical update: Mk2 modifies the lower diagonal
        % block of the next row before that row is used as a pivot.
        D{k+1}  = D{k+1}  - Mk2 * B{k};
        A0{k+2} = A0{k+2} - Mk2 * C{k};
        bb{k+2} = bb{k+2} - Mk2 * bb{k};
    end
end

% ===============================
% BACK SUBSTITUTION
% ===============================
x_cell = cell(nB,1);
x_cell{nB} = A0{nB} \ bb{nB};
x_cell{nB-1} = A0{nB-1} \ (bb{nB-1} - B{nB-1} * x_cell{nB});

for k = nB-2:-1:1
    x_cell{k} = A0{k} \ (bb{k} - B{k}*x_cell{k+1} - C{k}*x_cell{k+2});
end

x = vertcat(x_cell{:});

end