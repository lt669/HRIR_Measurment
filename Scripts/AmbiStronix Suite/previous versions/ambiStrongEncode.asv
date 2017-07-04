function varargout = ambistronixEncode(source, phi, theta, varargin)
    nChannels = nargout;
    varargout = cell(1,nChannels);
    
    norm = 'n3d'; % Default Normalisation
    
    for k = 1:length(varargin)
       switch varargin{k}
           case 'deg'
               phi = phi * pi / 180;
               theta = theta * pi / 180;
           case 'rad'
               % Do Nothing
           case 'n3d'
               norm = 'n3d';
           case 'sn3d'
               norm = 'sn3d';
           case 'orthonormal'
               norm = 'orthonormal';
           otherwise
               error('Invalid options')
       end
    end
    
    [L, M] = deal(0);
    
    for k = 1:nChannels;
        
        varargout{k} = ShHarmCoef_v1_0(L, M, phi, theta, norm) * source;

        M = M+1;
        
        if M > L
            L=L+1;
            M = -L;
        else
            %do nothing
        end
   end
    
end