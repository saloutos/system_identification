function Arm6Link_ellipsoid_visualize(obj, Phi, color_)
    
    % TODO: change joint angles here
    q = zeros(1,6);

    [model, ~]     = Arm6LinkModel();
    NB = model.NB;
    p  = model.parent;
    
    for i = 1 : NB
       I{i}    = inertiaVecToMat(Phi(:,i));
       Irot{i} = inertiaVecToMat(Phi(:,NB+i));
    end

    color = color_(1:NB,:);
    color_rot = color_(NB+1:end,:);
    
%     q_CAD = q + [0;pi;pi]; % Different joint angles to account for different zero config in CAD
    
    hold on;
    for i = 1:NB
      [ XJ, ~]     = jcalc( model.jtype{i}, q(i) );
      [ XJrot, ~]  = jcalc( model.jtype{i}, q(i)*model.gr{i} );
%       [ XJ_CAD, ~] = jcalc( model.jtype{i}, q_CAD(i) );
      
%       if i == 2
%          % Account for offset CAD origin for thigh link
%          XJ_CAD = XJ_CAD * [eye(3,3),zeros(3,3);skew([0;-0.045;0]),eye(3,3)]; 
%       end
      
      Xup{i} = XJ * model.Xtree{i};
      Xrot{i} = XJrot * model.Xrotor{i};
%       Xup_CAD{i} = XJ_CAD * model.Xtree{i};
      
      if model.parent(i) == 0
        X0{i}     = Xup{i};
        X0rot{i}  = Xrot{i};
%         X0_CAD{i} = Xup_CAD{i};
        
      else
        X0{i}     = Xup{i} *  X0{p(i)};
        X0rot{i}  = Xrot{i} * X0{p(i)};
%         X0_CAD{i} = Xup_CAD{i} *  X0{p(i)};
        
      end

      T    = inv( AdjointToSE3( X0{i} ) );
      Trot = inv( AdjointToSE3( X0rot{i} ) );

      plot_inertiatensor(T, I{i}, 0.45, color(i,:));
      plot_inertiatensor(Trot, Irot{i}, 0.45, color_rot(i,:));
        
      draw_SE3(T);
      
%       T = inv( AdjointToSE3( X0_CAD{i} ) );

      if i>0 %1
          V = (T(1:3,1:3)*obj{i}.v' + T(1:3,4) * ones(1,size(obj{i}.v,1)))';
          F = obj{i}.f.v;
          patch('vertices', V, 'faces', F,'LineStyle','none','FaceAlpha',.15);
      end
    end
    
    % add visualization for base link here too
    base_offset = [0,0,-0.051]';
    V = (eye(3)*obj{7}.v' + base_offset*ones(1,size(obj{7}.v,1)))';
    F = obj{7}.f.v;
    patch('vertices', V, 'faces', F,'LineStyle','none','FaceAlpha',.15);

    axis equal;
    view([-140 26]);
    draw_SE3(eye(4,4));
end