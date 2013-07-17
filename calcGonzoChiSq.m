function [chiSq] = calcGonzoChiSq(UP,GP,bias,N,var,...
    Gauss,GaussUni,Percents,Sigmas,varargin)

% This function is designed as part of the S2SconvoMahal.m function, but
% could be used for other purposes.  This function determines the
% appropriate Uniform Non-Central (Gonzo) Chi-Squared value:
%
% Inputs:
%   1.bias -> given a maximum bias
%   2.N    -> the number of input tracks to the convoluded gaussian
%   3.UP   -> the desired percentile of the convoluded gaussian to uniformly distribute over,
%   4.GP   -> the desired percentile of the Gonzo distribution.
%   5.var  -> the input covariance for bias vs. sigma calculations
%   6.Distribution Inputs
%     Name           Size                    Bytes  Class
%     __________________________________________________________
%     Gauss          1x25  -> Sigma points of a Gaussian corresponding to
%                             the Percents input
%     GaussUni      90x25  -> Sigma points of a Gaussian & Uniform
%                             distribution convolution correspoding to the
%                             Percents & Sigmas inputs
%     Percents       1x25  -> Percentile index [0:.5:99.5]
%     Sigmas         1x90  -> Sigma index [0:20]
%
%   7.varargin -> Region over which uncommon tracks are distributed
%
% Assumptions:
%   1. bias distribution
%   2. DoF
%   3. bias vs. convolution
%   4. file location
%
% Outputs:
%   1.chiSq -> DUH!
%   Note: chi square is an operation that acts on a variable, which has
%   traditionally been a Gaussian distribution.  However, now that
%   different distributions are convolved, this value is not static.
%
%
%==========================================================================

bias = (N*bias);
L = length( varargin);
% Here you are either setting the uncommon distribution against some known
% area, or to the standard Gaussian distribution
if L == 2
    a = varargin{1};
    b = varargin{2};
    Nuncomm = ceil(UP*N/100);
    % Using CLT, this is the sigma value of the uncommon track distribution
    UncommSigma = sqrt((b-a)^2/12 *Nuncomm);

    % Figure out the convoluded gaussian sigma; remember two track sets
    Csigma = sqrt(var*N*2);
else
    ix = find(Percents >=UP);
    Csigma = Gauss(ix(1));
end

% The bias is treated as the maximum possible bias between the two sensors,
% meaning that the bias could be anything inbetween.  So in effect, the
% maximum is added to the sigma value associated with the desired uniform
% percentile...except that the bias need to be in units of sigma, which
% depends specifically on the combined convoluded gaussians...

Bsigma = bias ./ sqrt(var(:));

% Combine the two sigmas
sigma = Csigma + Bsigma;
ixG = find(Gauss >=sigma);
if ~isempty( ixG)
    ixG = ixG(1); else
    ixG = length( Gauss); end
UPX2 = Gauss(ixG);

% OLD CODE OLD CODE
% ix = find(Gauss >=sigma);
% ix = ix(1);
% Pcnt = Percents(ix);
% Pcnt = UP + BP; % You cant add the percents together because they can
% easily go above 100, which doesnt make any sense, however the sigmas can
% go on to infinity (theoretically) and still be translated properly.
%
% OK, I AM TOO TIRED TO FIX THIS PROPERLY, BUT I SHOULD BE ABLE TO COMBINE
% THESE TWO INDEX LOOK UPS.  I ONLY DID THESE TWO STEPS BECAUSE I WAS
% ADDING THE TWO PERCENTS TOGETHER AND I SHOULDN'T BE.
%
% ixG = find(Percents >= Pcnt);
% if ~isempty( ixG)
%     ixG = ixG(1); else
%     ixG = length( Percents); end
% UPX2 = Gauss(ixG);
% OLD CODE OLD CODE

% The percentile is indexed with the associated sigma value, and since
% Gauss has only one hypothesis, central, it is only one dimension.
sV = length(var);
sS = length(Sigmas);
% UPX2 = Gsigma + sigma;

ixUP = find( UPX2 > max(Sigmas));
if ~isempty(ixUP)
    UPX2(ixUP) = max(Sigmas); end

SS = repmat(Sigmas(:),1,sV);
[ixSr,ixSc] = find( SS >= repmat(UPX2(:),sS,1));
xx = inf(sS,sV);
xx(ixSr,ixSc) = 1;
xx = repmat([1:sS],1,sV).*xx;
ixSr = min(xx,[],1); % Want the Sigma value immediately above desired threshold

ixP = find(Percents >= GP);
ixP = ixP(1);
if UPX2 >0
    chiSq = GaussUni(ixSr,ixP); else % Thresholds to use
    chiSq = Gauss(ixP)^2; end
