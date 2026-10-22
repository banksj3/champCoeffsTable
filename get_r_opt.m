function [rL,rR,kappa,A,info] = get_r_opt( theta,beta,lambda,method,file )
%% Optimal CHAMP interface coefficients r_L, r_R in [0,1] (order-2 / BDF2 scheme)
%% from the precomputed table.
%%
%%   [rL,rR]                  = get_r_opt( theta,beta,lambda )
%%   [rL,rR,kappa,A,info]     = get_r_opt( theta,beta,lambda,method,file )
%%
%% Inputs (arrays of equal size, or scalars):
%%   theta  = K_L/K_R,  beta = D_L/D_R,  lambda = dt*D_L/h^2
%%   method : 'linear' (default) or 'nearest'
%%   file   : table (default champ_order2.h5 next to this function)
%%
%% The interface condition on each side is (1-r)*N[T] + r*D[T] (r = 0 Neumann,
%% r = 1 Dirichlet); the Robin weight is p = r/(1-r) (see get_p_opt).
%%
%% The table holds several families of optima (local minimax solutions, e.g. A, A',
%% C, C'). Each family is evaluated separately, then the best one is returned:
%%   'linear'  : trilinear interpolation in log10(theta), log10(beta), log10(lambda)
%%               per family, of r in the measure log(r+0.01) and of the tracked and
%%               verified max|A|. A family takes part only if it is stable at all
%%               surrounding table nodes that carry weight (at a table node: that node). Among those, the family chosen at the
%%               nearest node is used if it takes part, otherwise the selection rule
%%               below. The interpolated r has not itself been verified; where no
%%               family takes part (e.g. next to the stability limit) the nearest-node
%%               value is returned instead (info.method = 'nearest' there).
%%   'nearest' : the nearest table node, with the family chosen there. Every value
%%               returned has been verified at that node.
%% Selection rule: smallest tracked minimax value A* unless several lie within
%% max(1.05 A*, A*+0.02) (or 0.1 times the verified max|A| when A* has collapsed);
%% among those the smallest verified max|A|, with ties resolved toward the family
%% chosen at neighbouring nodes.
%%
%% Outputs:
%%   rL, rR : coefficients; NaN where no family is stable or the query is outside
%%            the table (theta, beta in [1e-6,1e6], lambda in [1,1e8])
%%   kappa  : [size(theta) x 4] wavenumbers of the (up to 4) tracked peaks of
%%            |A(kappa)| of the chosen family, largest |A| first (NaN = no peak),
%%            from the table node nearest the query (an estimate for 'linear')
%%   A      : [size(theta) x 4] complex amplification factors of those peaks
%%   info   : Amax   verified max|A| over all roots and kappa (interpolated for 'linear')
%%            Astar  tracked minimax value A*
%%            family index of the chosen family (0 none stable, NaN outside)
%%            families  family names
%%            active [size x 4] peak is in the minimax active set
%%            method per query, 'linear' or 'nearest'

  persistent T
  if( nargin < 4 || isempty(method) ), method = 'linear'; end
  if( nargin < 5 || isempty(file) ), file = fullfile( fileparts(mfilename('fullpath')),'champ_order2.h5' ); end
  if( isempty(T) || ~strcmp(T.file,file) ), T = loadTable( file ); end
  sz = size(theta); n = numel(theta);
  q = [log10(theta(:)) log10(beta(:)) log10(lambda(:))].';
  g = {T.lt,T.lb,T.ll}; nF = numel(T.names);
  inR = true(1,n);
  for d = 1:3, inR = inR & isfinite(q(d,:)) & q(d,:) >= g{d}(1)-1e-12 & q(d,:) <= g{d}(end)+1e-12; end

  % nearest node (always needed: family choice, peaks, fallback)
  idx = ones(3,n);
  for d = 1:3, idx(d,:) = min( max( round((q(d,:)-g{d}(1))/(g{d}(2)-g{d}(1)))+1,1 ),numel(g{d}) ); end
  nn = sub2ind( T.sz,idx(1,:),idx(2,:),idx(3,:) );
  bestN = T.best(nn); bestN(~inR) = 0;

  useLin = strcmpi(method,'linear');
  if( ~useLin && ~strcmpi(method,'nearest') ), error( 'get_r_opt: unknown method %s',method ); end
  fam = bestN; rLq = nan(1,n); rRq = rLq; Am = rLq; As = rLq; how = repmat( {'nearest'},1,n );
  if( useLin )
    lo = ones(3,n); w = zeros(3,n);
    for d = 1:3
      s = (q(d,:)-g{d}(1))/(g{d}(2)-g{d}(1));
      lo(d,:) = min( max( floor(s)+1,1 ),numel(g{d})-1 ); w(d,:) = min( max( s-(lo(d,:)-1),0 ),1 );
    end
    bl = @(r) log(r+0.01);
    RL = nan(nF,n); RR = RL; AT = RL; AM = RL; ok = false(nF,n);
    for f = 1:nF
      B = T.F(f); acc = zeros(4,n); okf = true(1,n);
      for c = 0:7
        b = bitget( c,1:3 ).';
        li = sub2ind( T.sz,lo(1,:)+b(1),lo(2,:)+b(2),lo(3,:)+b(3) );
        wc = prod( b.*w+(1-b).*(1-w),1 );
        on = wc >= 1e-12;                                  % corners with zero weight do not count
        okf = okf & ( B.stable(li) | ~on );
        v = [bl(B.rL(li)); bl(B.rR(li)); B.Atracked(li); B.Amax(li)]; v(:,~on) = 0;
        acc = acc+wc.*v;
      end
      RL(f,:) = min( max(exp(acc(1,:))-0.01,0),1 ); RR(f,:) = min( max(exp(acc(2,:))-0.01,0),1 );
      AT(f,:) = acc(3,:); AM(f,:) = acc(4,:); ok(f,:) = okf & inR;
    end
    sel = selectFamily( AT,AM,ok );
    useN = bestN > 0; useN(useN) = ok( sub2ind(size(ok),bestN(useN),find(useN)) );   % the nearest node's family, where valid
    sel(useN) = bestN(useN);
    j = find( sel > 0 ); li = sub2ind( size(RL),sel(j),j );
    rLq(j) = RL(li); rRq(j) = RR(li); Am(j) = AM(li); As(j) = AT(li); fam(j) = sel(j); how(j) = {'linear'};
  end
  % nearest-node values where no interpolated family was available
  j = find( isnan(rLq) & bestN > 0 );
  for i = j
    B = T.F(bestN(i)); rLq(i) = B.rL(nn(i)); rRq(i) = B.rR(nn(i)); Am(i) = B.Amax(nn(i)); As(i) = B.Atracked(nn(i)); fam(i) = bestN(i);
  end
  fam(isnan(rLq)) = 0;

  % peaks of the chosen family at the nearest node
  K = nan(n,4); Ap = complex(K,K); act = false(n,4);
  for i = find( fam > 0 )
    B = T.F(fam(i));
    K(i,:) = B.peakKappa(:,nn(i)).'; Ap(i,:) = complex( B.peakAreal(:,nn(i)),B.peakAimag(:,nn(i)) ).'; act(i,:) = B.peakActive(:,nn(i)).' > 0;
  end

  fam = double(fam); fam(~inR) = NaN;
  if( any(~inR) ), warning( 'get_r_opt:range','%d of %d queries outside the table range: NaN returned',sum(~inR),n ); end
  rL = reshape( rLq,sz ); rR = reshape( rRq,sz );
  kappa = reshape( K,[sz 4] ); A = reshape( Ap,[sz 4] );
  if( isvector(theta) || isscalar(theta) ), kappa = reshape( K,[n 4] ); A = reshape( Ap,[n 4] ); end
  info.Amax = reshape( Am,sz ); info.Astar = reshape( As,sz ); info.family = reshape( fam,sz );
  info.families = T.names; info.active = act; info.method = reshape( how,sz );
end

function best = selectFamily( Atr,Amax,ok )
  % rows: families, columns: queries
  Atr(~ok) = Inf;
  [Amin,q] = min( Atr,[],1 ); cap = max( 1.05*Amin,Amin+0.02 );
  j = find( isfinite(Amin) ); Am0 = Amax( sub2ind(size(Amax),q(j),j) );
  coll = Amin(j) < 0.1*Am0; cap(j(coll)) = max( cap(j(coll)),0.1*Am0(coll) );
  M = Amax; M( ~(ok & Atr <= cap) ) = Inf;
  [m,best] = min( M,[],1 ); best( ~isfinite(m) ) = 0;
end

function T = loadTable( file )
  T.file = file;
  T.lt = log10( h5read(file,'/theta') ); T.lb = log10( h5read(file,'/beta') ); T.ll = log10( h5read(file,'/lambda') );
  T.lt = T.lt(:).'; T.lb = T.lb(:).'; T.ll = T.ll(:).';
  T.sz = [numel(T.lt) numel(T.lb) numel(T.ll)];
  T.best = h5read( file,'/best' );
  T.names = strsplit( h5readatt(file,'/','families'),',' );
  for f = 1:numel(T.names)
    for v = {'rL','rR','Atracked','Amax','status','peakKappa','peakAreal','peakAimag','peakActive'}
      T.F(f).(v{1}) = h5read( file,['/' T.names{f} '/' v{1}] );
    end
    T.F(f).stable = T.F(f).status > 0 & T.F(f).Amax <= 1+1e-6;
    T.F(f).peakKappa = reshape( T.F(f).peakKappa,4,[] ); T.F(f).peakAreal = reshape( T.F(f).peakAreal,4,[] );
    T.F(f).peakAimag = reshape( T.F(f).peakAimag,4,[] ); T.F(f).peakActive = reshape( T.F(f).peakActive,4,[] );
  end
end
