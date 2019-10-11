classdef BWD
 
 properties(Constant)
  sig_ = 2 .^ (2:9);
  tau_ = 2 .^ (1:2);
  szf  = 9;%WARNING very sensitive!!!!!
  tszf = 9;
  fields = BWD.getRegions();
  o = cos([0,pi/2,pi,3*pi/2]) + 1j * sin([0,pi/2,pi,3*pi/2]);
 end
 
 methods(Static = true)
  
  %% Binary String generation
  
  function code = encode(v,stip)
   
   [sv,s_,t_] = BWD.subvid(v,stip);
   
   scan = BWD.fields{s_,t_};
   
   sv = BWD.flipvideo(sv);
   
   sv = double(sv(:));
   
   code = mexbwd(sv,scan{1},scan{2});
   
  end
  
  function [sv,s_,t_] = subvid(v,stip)
   stip = double(stip);
   s_ = floor(BWD.szf  * sqrt(stip(4)));
   t_ = floor(BWD.tszf * sqrt(stip(5)));
   
   y0 = stip(1);
   x0 = stip(2);
   t0 = stip(3);
   
   [Y,X,T] = size(v);
   
   ry = round((1:s_)+y0-s_/2);
   rx = round((1:s_)+x0-s_/2);
   rt = round((1:t_)+t0-t_/2);
   
   ry(ry < 1) = ry(1);
   ry(ry > Y) = ry(end);
   
   rx(rx < 1) = rx(1);
   rx(rx > X) = rx(end);
   
   rt(rt < 1) = rt(1);
   rt(rt > T) = rt(end);
   
   sv = v(ry,rx,rt);
  end
  
  function vr = flipvideo(v)
   vr = v;
   v = abs(v);
   [Y,X,T] = size(v);
   hv = v(:,:,1:floor(T/2));
   [r,~] = BWD.xytSpV(Y,X,1);
   hv = sum(hv,3);
   hv = hv(:);
   x = hv .* r(:,1);
   y = hv .* r(:,2);
   c = sum([y,x])/(sum(hv)) - [Y,X]/2;
   c = c(1)+1j*c(2);
   c = c/abs(c);
   c = BWD.o-c;
   c = c .* conj(c);
   [~,c] = min(c);
   switch c
    case 1
     vr =  fliplr(vr);
    case 2
     vr =  flipud(vr);
    otherwise
    return
   end
  end
  
  %% Scanning Regions
  
  function v = fillvolumeidx(v)
   idx = 0;
   for t = 1:size(v,3)
    for x = 1:size(v,2)
     for y = 1:size(v,1)
      idx = idx+1;
      sv = v{y,x,t};
      sv(:,:,:) = idx;
      v{y,x,t} = sv;
     end
    end
   end
   v = cell2mat(v);
   v = v(:);
  end
    
  function regions = getRegions(varargin)
   
   if nargin == 1
    st  = varargin{1};
    ts_ = st(1);
    ss_ = st(2);
   else
    ts_ = BWD.tau_;
    ss_ = BWD.sig_;
   end
   regions = cell(2,2);
   
   for t = ts_
    for s = ss_
     s_ = floor(BWD.szf  * sqrt(s));
     t_ = floor(BWD.tszf * sqrt(t));
     v = zeros(s_,s_,t_);
     
     Y = BWD.splitInt(size(v,1),4);
     X = BWD.splitInt(size(v,2),4);
     T = BWD.splitInt(size(v,3),2);
     
     v = mat2cell(v,Y,X,T);
     
     v = BWD.fillvolumeidx(v);
     
     pttn = BWD.pttnfield();
     scan = cell(1e2,2);
     
     n = 0;
     
     %space symmetric
     for k = 1:32%size(pttn,1)
      A = [pttn{k,1};pttn{k,1}+16];
      B = [pttn{k,2};pttn{k,2}+16];
      n = n+1;
      x = ismember(v,A);
      y = ismember(v,B);
      scan(n,:) = {x,y};
     end
     
     % time symmetric
     %%{
     for k = 1:32%size(pttn,1)
      A = [pttn{k,1};pttn{k,2}+16];
      B = [pttn{k,2};pttn{k,1}+16];
      n = n+1;
      x = ismember(v,A);
      y = ismember(v,B);
      scan(n,:) = {x,y};
     end
       
     scan = scan(1:n,:)';
     
     scan = {cell2mat(scan(1,1:64)),cell2mat(scan(2,1:64))};
    
     regions{s_,t_} = scan;
     
    end
   end
   
   if numel(ts_) == 1
    regions = {scan};
   end
  
  end
  
  function [r,p] = xytSpV(nx,ny,nt)
   p = struct;
   p.x = zeros(nx,ny,nt);
   p.y = zeros(nx,ny,nt);
   p.t = zeros(nx,ny,nt);
   for t = 1:nt
    for x = 1:nx
     for y = 1:ny
      p.y(y,x,t) = y;
      p.x(y,x,t) = x;
      p.t(y,x,t) = t;
     end
    end
   end
   a = p.y;
   b = p.x;
   c = p.t;
   r = [a(:),b(:),c(:)];
  end
 
  function field = pttnfield()
   field = cell(1,2);
   field(01,:) = {[1;2;3;4;5;6;7;8],setdiff((1:16)',[1;2;3;4;5;6;7;8])};
   field(02,:) = {[1;5;9;13;2;6;10;14],setdiff((1:16)',[1;5;9;13;2;6;10;14])};
   field(03,:) = {[1;5;9;13;4;8;12;16],setdiff((1:16)',[1;5;9;13;4;8;12;16])};
   field(04,:) = {[1;2;3;4;13;14;15;16],setdiff((1:16)',[1;2;3;4;13;14;15;16])};
   field(05,:) = {[1;2;3;4;9;10;11;12],setdiff((1:16)',[1;2;3;4;9;10;11;12])};
   field(06,:) = {[1;5;9;13;3;7;11;15],setdiff((1:16)',[1;5;9;13;3;7;11;15])};
   field(07,:) = {[1;2;7;8;11;12;13;14],setdiff((1:16)',[1;2;7;8;11;12;13;14])};
   field(08,:) = {[2;6;3;7;9;13;12;16],setdiff((1:16)',[2;6;3;7;9;13;12;16])};
   field(09,:) = {[1;5;2;6;11;15;12;16],setdiff((1:16)',[1;5;2;6;11;15;12;16])};
   field(10,:) = {[1;4;6;7;10;11;13;16],setdiff((1:16)',[1;4;6;7;10;11;13;16])};
   % sparse 8
   %%{
   field(11,:) = {[1;5;3;7;10;14;12;16],setdiff((1:16)',[1;5;3;7;10;14;12;16])};%corrected
   field(12,:) = {[1;2;7;8;9;10;15;16],setdiff((1:16)',[1;2;7;8;9;10;15;16])};
   field(13,:) = {[1;3;6;10;8;12;13;15],setdiff((1:16)',[1;3;6;10;8;12;13;15])};%corrected
   field(14,:) = {[1;4;6;7;9;12;14;15],setdiff((1:16)',[1;4;6;7;9;12;14;15])};%corrected
   
   field(15,:) = {[1;3;6;8;9;11;14;16],setdiff((1:16)',[1;3;6;8;9;11;14;16])};
   %}
   
   field(16,:) = {[1;2;5;6],[3;7;4;8]};
   field(17,:) = {[9;13;10;14],[11;15;12;16]};
   field(18,:) = {[3;7;4;8],[11;15;12;16]};
   field(19,:) = {[1;5;2;6],[9;13;10;14]};
   
   field(20,:) = {[1;2;3;4],[5;6;7;8]};
   field(21,:) = {[3;7;11;15],[4;8;12;16]};
   field(22,:) = {[9;10;11;12],[13;14;15;16]};
   field(23,:) = {[1;5;9;13],[2;6;10;14]};
   
   field(24,:) = {[1;5;3;7],[2;6;4;8]};
   field(25,:) = {[3;4;11;12],[7;8;15;16]};
   field(26,:) = {[9;13;11;15],[10;14;12;16]};
   field(27,:) = {[1;2;9;10],[5;6;13;14]};
   
   field(28,:) = {[1;5;2;6],[11;15;12;16]};
   field(29,:) = {[3;7;4;8],[9;10;13;14]};
   
   field(30,:) = {[2;3;14;15],[5;9;8;12]};
   
   field(31,:) = {[1;5;4;8],[10;11;14;15]};
   field(32,:) = {[1;2;13;14],[7;11;8;12]};
   %89.35% g 1 c 500
   
   %{
   field(33,:) = {[1;4;13;16],[6;7;10;11]};
   field(34,:) = {[2;8;9;15],[3;5;12;14]};
   
   field(35,:) = {[2;6;3;7],[9;13;12;16]};
   field(36,:) = {[5;6;9;10],[3;4;15;16]};
   
   field(37,:) = {[2;3;6;7],[10;11;14;15]};
   field(38,:) = {[5;6;9;10],[7;8;11;12]};
   field(39,:) = {[2;3;14;15],[6;7;10;11]};
   field(40,:) = {[5;9;8;12],[6;7;10;11]};

   field(41,:) = {[1;2;3;4],[13;14;15;16]};
   field(42,:) = {[1;5;9;13],[4;8;12;16]};
   field(43,:) = {[1;5;12;16],[4;8;9;13]};
   field(44,:) = {[1;2;15;16],[3;4;13;14]};

   field(45,:) = {[3;8;12;15],[4;7;11;16]};
   field(46,:) = {[9;12;14;15],[10;11;13;16]};
   field(47,:) = {[2;5;9;14],[1;6;10;13]};
   field(48,:) = {[2;3;5;8],[1;4;6;7]};

   field(49,:) = {[2;6;10;14],[3;7;11;15]};
   field(50,:) = {[5;6;7;8],[9;10;11;12]};
   field(51,:) = {[2;6;11;15],[3;7;10;14]};
   field(52,:) = {[5;6;11;14],[7;8;9;10]};

   field(53,:) = {[1;5;11;15],[3;7;9;13]};
   field(54,:) = {[1;2;11;12],[3;4;9;10]};
   field(55,:) = {[2;6;12;16],[4;8;10;14]};
   field(56,:) = {[5;6;15;16],[7;8;14;14]};

   field(57,:) = {[2;7;11;14],[3;6;10;15]};
   field(58,:) = {[5;8;10;11],[6;7;9;12]};
   field(59,:) = {[5;8;9;12],[1;4;13;16]};
   field(60,:) = {[2;3;14;15],[1;4;13;16]};

   field(61,:) = {[1;5;3;7],[9;13;11;15]};
   field(62,:) = {[1;2;9;10],[3;4;11;12]};
   field(63,:) = {[2;4;6;8],[10;12;14;16]};
   field(64,:) = {[5;6;13;14],[7;8;15;16]};
   %}
   
   
   %{
   %26 X4
   
   field(100,:) = {[2;5;10;13],[1;6;9;14]};
   field(101,:) = {[1;6;3;8],[2;5;7;4]};
   field(102,:) = {[4;7;12;15],[3;8;11;16]};
   field(103,:) = {[9;14;11;16],[10;12;13;15]};
   
   field(104,:) = {[1;4;6;7],[2;3;5;8]};
   field(105,:) = {[3;8;12;15],[4;8;11;16]};
   field(106,:) = {[9;12;14;15],[10;11;13;16]};
   field(107,:) = {[1;3;5;7],[9;11;13;15]};
   
   field(108,:) = {[1;2;9;10],[3;4;11;12]};
   field(109,:) = {[2;6;4;8],[10;12;14;16]};
   field(110,:) = {[5;6;13;14],[7;8;15;16]};
   field(111,:) = {[2;7;11;14],[3;6;10;15]};
   
   field(112,:) = {[5;8;10;11],[6;7;9;12]};
   field(113,:) = {[1;6;10;13],[2;5;9;14]};
   field(114,:) = {[6;7;10;11],[5;8;9;12]};
   field(115,:) = {[2;3;14;15],[6;7;10;11]};
   
   field(116,:) = {[2;8;9;15],[3;5;12;14]};
   field(117,:) = {[1;4;13;16],[6;7;10;11]};
   field(118,:) = {[1;4;13;16],[2;3;14;15]};
   field(119,:) = {[1;4;13;16],[5;8;9;12]};
   
   field(120,:) = {[1;4;9;12],[5;8;13;16]};
   field(121,:) = {[1;3;13;15],[2;4;14;16]};
   field(122,:) = {[3;4;15;16],[5;6;9;10]};
   field(123,:) = {[2;3;6;7],[9;12;13;16]};
   
   field(124,:) = {[5;7;10;12],[6;8;9;11]};
   field(125,:) = {[2;7;10;15],[3;6;11;14]};

   % time symmetric
   field(126,:) = {[1;2;5;6]};
   field(127,:) = {[3;4;7;8]};
   field(128,:) = {[9;10;13;14]};
   field(129,:) = {[11,12,15,16]};
   field(130,:) = {[9;10;11;12;13;14;15;16]};
   field(131,:) = {[1;2;5;6;9;10;13;14]};
   
   field(132,:) = {[1;2;3;4;5;6;7;8]};
   field(133,:) = {[3;4;7;8;11;12;15;16]};
   field(134,:) = {[1;4;13;16]};
   field(135,:) = {[6;7;10;11]};
   field(136,:) = {(1:16)'};
   field(137,:) = {(1:16)'};
   %}
   

   
   for k = size(field,1):-1:1
    if isempty(field{k,1})
     field(k,:) = [];
    end
   end
   
  end
  
  function s = splitInt(a,b)
   x = floor(a/b);
   s = x * ones(b,1);
   m = 0; 
   
   while sum(s) ~= a
     m = m+1;
     s(m) = s(m)+1;
   end
  end
  
 end
end