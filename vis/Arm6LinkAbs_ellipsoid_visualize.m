function Arm6LinkAbs_ellipsoid_visualize(obj, Phi, color_)
    
    % TODO: change joint angles here
    q = zeros(1,6);

    [model, ~]     = Arm6LinkModel();
    N_RB = model.N_RB;
    NJ = model.NB;
    NL = model.NL;
    p  = model.parent;
    
    for i = 1 : NL % need to make sure Phi is in the correct order?
       I{i}    = inertiaVecToMat(Phi(:,i));
       Irot{i} = inertiaVecToMat(Phi(:,NL+i));
    end

    color = color_(1:NL,:);
    color_rot = color_(NL+1:end,:);
        
    hold on;
    for i = 1:NL % for each joint

      % TODO: probably not handling absolute pair joint correctly
      if i==4 % if i is 4 (first link of abs pair)
          [ XJ, ~]     = jcalc( ['R' model.joint{i}.jointAxis{1}], q(i) );
          n1 = model.joint{i}.gearRatio{1}*model.joint{i}.beltRatio{1};
          [ XJrot, ~]  = jcalc( ['R' model.joint{i}.rotorAxis{1}], q(i)*n1);
      elseif i==5 % else if i is 5 (second link of abs pair)
          k = 4;
          [ XJ, ~]     = jcalc( ['R' model.joint{k}.jointAxis{2}], q(i) );
          b1 = model.joint{k}.beltRatio{1};
          gr2 = model.joint{k}.gearRatio{2};
          n2 = gr2*model.joint{k}.beltRatio{2};
          [ XJrot, ~]  = jcalc( ['R' model.joint{k}.rotorAxis{2}], n2*q(i) + gr2*b1*q(k));
      elseif i==6 % last link is actually 5th joint model
          k = 5;
          [ XJ, ~]     = jcalc( ['R' model.joint{k}.jointAxis], q(i) );
          [ XJrot, ~]  = jcalc( ['R' model.joint{k}.rotorAxis], q(i)*model.joint{k}.gearRatio);
      else % else (i is 1,2,3)
          [ XJ, ~]     = jcalc( ['R' model.joint{i}.jointAxis], q(i) );
          [ XJrot, ~]  = jcalc( ['R' model.joint{i}.rotorAxis], q(i)*model.joint{i}.gearRatio);
      end
         
      Xup{i} = XJ * model.XtreeKin{i};
      Xrot{i} = XJrot * model.XtreeKinRot{i};
      
      if i==1 %model.parent(i) == 0
        X0{i}     = Xup{i};
        X0rot{i}  = Xrot{i};    
      else
        X0{i}     = Xup{i} *  X0{i-1}; %X0{p(i)};
        if i==5
            X0rot{i}  = Xrot{i} * X0{i-2}; %X0{p(i)}; % parent link is one higher
        else
            X0rot{i}  = Xrot{i} * X0{i-1}; %X0{p(i)};
        end
      end

      T    = inv( AdjointToSE3( X0{i} ) );
      Trot = inv( AdjointToSE3( X0rot{i} ) );

      plot_inertiatensor(T, I{i}, 0.45, color(i,:));
      plot_inertiatensor(Trot, Irot{i}, 0.45, color_rot(i,:));
        
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