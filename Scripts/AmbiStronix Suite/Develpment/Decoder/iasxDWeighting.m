function Dw = iasxDWeighting( D, weighting, M )

% Version 1.0 Feb 2017
%
% [ Weighted Decoding Matrix ] From f( Decoding Matrix,
%                                      Weighting,
%                                      Degree )
%
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 14 / 03 / 2017: Initial coding

%% AMP. / ENERGY PRESERVATION NORMALISATION

    a_0 = 1; % Amp. Pres.
%   Energy Pres. ?

%%  DECLARATIONS
    syms fP_M(u) % Declare Legendre Polynomial (LP)
    
%% CODE
    for m = 0:M
       
        switch weighting
            case 'basic'
                a_m_Dash = 1;
                
            case 'maxre' % 3d
                fP_M(u) = (1 / (2^(M+1) * factorial((M+1))))*diff((u^2 - 1)^(M+1), u , (M+1)); %#ok<AGROW>
                roots = root(fP_M, u);
                maxrE = max(roots); 
                fP_M(u) = (1 / (2^m * factorial(m)))*diff((u^2 - 1)^m, u , m); %#ok<AGROW>
                a_m_Dash = fP_M(maxrE);    
                
            case 'inphase' % 3d
                a_m_Dash = (factorial(M) * factorial(M + 1)) / (factorial(M + m + 1) * factorial(M - m));
        end
        
        a_M_Dash((m)^2 + 1 : (m)^2 + 1 + (2*m)) =  double(a_m_Dash); %3d
        
    end
    
    a_M = a_M_Dash * a_0;
    Gamma = diag(a_M);
    Dw = D * Gamma;
    
end