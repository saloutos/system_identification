function Arm6Link_bound_visualize(obj, Q, color_)
    
    % TODO: change joint angles here
    q = zeros(1,6);

    [model, ~]     = Arm6LinkModel();
    N_RB = model.N_RB;
    NJ = model.NB;
    p  = model.parent;
    
    color = color_(1:NJ,:);
    color_rot = color_(NJ+1:end,:);
    
    hold on;
    for i = 1:NJ % for each joint
      [ XJ, ~]     = jcalc( ['R' model.joint{i}.jointAxis], q(i) );
      [ XJrot, ~]  = jcalc( ['R' model.joint{i}.rotorAxis], q(i)*model.joint{i}.gearRatio);
      
      Xup{i} = XJ * model.Xtree{i}(1:6,1:6);
      Xrot{i} = XJrot * model.Xtree{i}(7:12,1:6);
      
      if model.parent(i) == 0
        X0{i}     = Xup{i};
        X0rot{i}  = Xrot{i};    
      else
        X0{i}     = Xup{i} *  X0{p(i)};
        X0rot{i}  = Xrot{i} * X0{p(i)};
      end

      T    = inv( AdjointToSE3( X0{i} ) );
      Trot = inv( AdjointToSE3( X0rot{i} ) );

      plot_bounding_ellipsoid(T, Q{i}, .45, color(i,:))
      plot_bounding_ellipsoid(Trot, Q{i+NJ}, .45, color_rot(i,:))
      
      draw_SE3(T);

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

