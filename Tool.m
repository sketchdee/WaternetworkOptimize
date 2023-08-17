classdef Tool
    %TOOL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        Property1
    end
    
    properties(Constant)
        Property2=1;
    end
    
    %static methods
    methods(Static)
        
        %%
        %M2V用于将矩阵转换为行向向量，直接拼接
        %V2M用于将矩阵还原为向量， 也可以做
        function [M] = V2M(V,m,n)
            M=[];
            for i=1:1:m
                M=[M;V(1,(1+n*(i-1)):n*i)];
            end
        end
        
        function [V,m,n] = M2V(M) 
            [m,n]=size(M);
            V=[];
            for i=1:1:m
                V=[V,M(i,:)];
            end
        end
    end
    
    % common methods
    methods
        
    end
    
end

