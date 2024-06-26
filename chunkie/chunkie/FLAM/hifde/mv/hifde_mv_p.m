% HIFDE_MV_P   Dispatch for HIFDE_MV with F.SYMM = 'P'.

function X = hifde_mv_p(F,X)

  % initialize
  n = F.lvp(end);

  % upward sweep
  for i = 1:n
    sk = F.factors(i).sk;
    rd = F.factors(i).rd;
    T = F.factors(i).T;
    if ~isempty(T), X(sk,:) = X(sk,:) + T*X(rd,:); end
    X(rd,:) = F.factors(i).L'*X(rd,:);
    X(rd,:) = X(rd,:) + F.factors(i).E'*X(sk,:);
  end

  % downward sweep
  for i = n:-1:1
    sk = F.factors(i).sk;
    rd = F.factors(i).rd;
    T = F.factors(i).T;
    X(sk,:) = X(sk,:) + F.factors(i).E*X(rd,:);
    X(rd,:) = F.factors(i).L*X(rd,:);
    if ~isempty(T), X(rd,:) = X(rd,:) + T'*X(sk,:); end
  end
end