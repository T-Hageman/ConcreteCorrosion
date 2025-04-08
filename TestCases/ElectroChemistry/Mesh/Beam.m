clear all;
%close all;
clc
rng('default');

%This requires git clone https://github.com/YingzhouLi/meshpart.git
addpath(genpath('../../PostProcessingReqs/meshpart-master'))

	nCores = 32;
	fname = "Beam_32c";
	
	figure;
	t = tiledlayout('flow');
	t.Padding = "tight";
	t.TileSpacing = "tight";
	
	Lx = 0.1; 
	Ly = 0.05;
	Lz = 0.05;
	
	RBar = 5e-3;
	RPit = 4e-4;
	APit = 2*pi*RPit^2
	
	dxPit = 3*RPit/5;
	dxBar = 2*RBar/2.5;
	dxMax = 2*RBar*2;
	
	
	
	model = createpde(1);
	importGeometry(model,"Geo.stl");
	
	nexttile
	pdegplot(model,"FaceAlpha",0.3,'FaceLabels','off','CellLabels','off','VertexLabels','On')
	
	model = generateMesh(model,'Hgrad',1.2,"Hmax",dxMax,"Hface",{7,dxPit,6,dxBar});
	nexttile
	pdeplot3D(model,'NodeLabels','off','ElementLabels','off')
	disp("Model Volume: "+string(model.volume))
	
	% Nodes
		Nodes = model.Nodes';
		disp("Number of Nodes: "+string(length(Nodes)))
	
	%% interior elements
		Elementgroups{1}.Name = "Interior";
		Elementgroups{1}.Type = "T3_10B";
		Elementgroups{1}.Elements = model.Elements(:,findElements(model,'region','Cell',1))';
	
	nexttile
	xyz = Nodes(Elementgroups{1}.Elements(1,:),:);
	for n=1:10
		plot3(xyz(n,1),xyz(n,2),xyz(n,3),'k*')
		hold on
		text(xyz(n,1),xyz(n,2),xyz(n,3),"N"+string(n))
	end
	sets = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4];
	for i=1:6
		set = sets(i,:);
		plot3(xyz(set,1),xyz(set,2),xyz(set,3),'k')
	end
	
	%% Metal surface
		Elementgroups{2}.Name = "Pit";
		Elementgroups{2}.Type = "T2_6B";
		Elementgroups = getBoundaryGroup(Elementgroups, 2, model, 7);
	
		Elementgroups{3}.Name = "Bar";
		Elementgroups{3}.Type = "T2_6B";
		Elementgroups = getBoundaryGroup(Elementgroups, 3, model, 6);
	
	%% Exterior
		Elementgroups{4}.Name = "Top";
		Elementgroups{4}.Type = "T2_6B";
		Elementgroups = getBoundaryGroup(Elementgroups, 4, model, 4);
	
		Elementgroups{5}.Name = "Left";
		Elementgroups{5}.Type = "T2_6B";
		Elementgroups = getBoundaryGroup(Elementgroups, 5, model, 1);
	
	%% symmetry
 		Elementgroups{6}.Name = "Bottom";
		Elementgroups{6}.Type = "T2_6B";
		Elementgroups = getBoundaryGroup(Elementgroups, 6, model, 3);
	
		Elementgroups{7}.Name = "Right";
		Elementgroups{7}.Type = "T2_6B";
		Elementgroups = getBoundaryGroup(Elementgroups, 7, model, 8);
	
	%% Ends
		Elementgroups{8}.Name = "Symmetry";
		Elementgroups{8}.Type = "T2_6B";
		Elementgroups = getBoundaryGroup(Elementgroups, 8, model, 2);
	
		Elementgroups{9}.Name = "Back";
		Elementgroups{9}.Type = "T2_6B";
		Elementgroups = getBoundaryGroup(Elementgroups, 9, model, 5);
	
	nexttile
	PlotMesh(Nodes,Elementgroups,[2 3]);
	xlim([0.005 0.015])
	ylim([-0.005 0.005])
	zlim([-0.01 0])
	
	%% NodeGroups
	for eg=1:length(Elementgroups)
		Nodegroups{eg}.Name = Elementgroups{eg}.Name;
		Nodegroups{eg}.Nodes = sort(unique(Elementgroups{eg}.Elements(:)));
	end
	
	Nodegroups{eg+1}.Name = "Vertex";
	Nodegroups{eg+1}.Nodes = model.findNodes("nearest",[0,0,0]');
	
	%% Partitioning
	[Nodes, Nodegroups, Elementgroups, coreData] = partition(Nodes, Nodegroups, Elementgroups, nCores);
	
	nexttile
	PlotPartition(Nodes,Elementgroups,coreData);
	
	drawnow();
	saveToHDF(fname, Nodes, Nodegroups, Elementgroups, coreData)



function PlotPartition(nodes, ElementGroups, coreData)
	clrs = distinguishable_colors(coreData.NCores);
	if true
		for f=1:4
			switch f
				case 1
					order = [1 2 3];
				case 2
					order = [1 2 4];
				case 3
					order = [2 3 4];
				case 4
					order = [3 1 4];
			end
			F = ElementGroups{1}.Elements(:,order);
			CNum = coreData.ElementGroup{1}.core(:);
			C = clrs(CNum,:);
			%C(CNum~=7,:)=NaN;
			patch('Faces',F,'Vertices',nodes,'FaceVertexCData',C,'EdgeColor','k','FaceColor','flat');
			hold on
		end
	else
		for eg=2:length(ElementGroups)
			if ElementGroups{eg}.Type == "T2_6"
				%for el=1:length(ElementGroups{eg}.Elements)
					F = ElementGroups{eg}.Elements(:,[1 2 3]);
					C = clrs(coreData.ElementGroup{eg}.core(:),:);
					%patch('Faces',F,'Vertices',nodes,'FaceColor',C,'FaceAlpha',.5)
					patch('Faces',F,'Vertices',nodes,'FaceVertexCData',C,'EdgeColor','k','FaceColor','flat');
					hold on
				%end
			end
		end
	end
end

function [Nodes, NodeGroup, ElementGroup, coreData] = partition(Nodes, NodeGroup, ElementGroup, nCores)
	coreData.NCores = nCores;

	NNodes = length(Nodes);

	ix = [];
	iy = [];
	iz = [];
	for el=1:length(ElementGroup{1}.Elements)
		nds = ElementGroup{1}.Elements(el,:);
		[allocX, allocY] = ndgrid(nds,nds);
		allocZ =  allocX*0+1;
		ix = [ix allocX(:)];	
		iy = [iy allocY(:)];	
		iz = [iz allocZ(:)];	
	end
	Connections = sparse(ix, iy, iz, NNodes, NNodes);

	map = geodice(Connections,Nodes,log2(nCores));
	%map = specdice(Connections,log2(nCores));
	coreData.NodeLocs = map'+1;

	%NodeGroups
	for eg=1:length(NodeGroup)
		for el=1:length(NodeGroup{eg}.Nodes)
			coreData.NodeGroup{eg}.core(el) = coreData.NodeLocs(NodeGroup{eg}.Nodes(el));
		end
	end

	BeenClaimed = zeros([length(Nodes),1]);

	%Elements
	for eg=1:length(ElementGroup)
		for el=1:size(ElementGroup{eg}.Elements,1)
			eNodes = ElementGroup{eg}.Elements(el,:);
			cNodes = coreData.NodeLocs(eNodes);
			loc = mode(cNodes);
			coreData.ElementGroup{eg}.core(el) = loc;
		end
	end

	% private/shared Nodes
	for c=1:nCores
		coreData.Nodes.PrivateNodes{c} = find(coreData.NodeLocs == c);

		coreData.Nodes.Ghosts_ImReceiving{c} = [];
		coreData.Nodes.Ghosts_ImSending{c} = [];
	end
	
	for eg=1:length(ElementGroup)
		for el=1:size(ElementGroup{eg}.Elements,1)
			Nds = ElementGroup{eg}.Elements(el,:);
			NodeCores = coreData.NodeLocs(Nds);
			uniqCores = sort(unique(NodeCores));
			if (length(uniqCores)>1)
				ElemOwner = coreData.ElementGroup{eg}.core(el);
				Ghosts = find(NodeCores~=ElemOwner);
				GhostNodes = Nds(Ghosts);
				GhostCores = NodeCores(Ghosts);
				for j=1:length(Ghosts)
					coreData.Nodes.Ghosts_ImReceiving{ElemOwner}(end+1) = GhostNodes(j);
					coreData.Nodes.Ghosts_ImSending{GhostCores(j)}(end+1) = GhostNodes(j);
				end
			end
		end
	end
	for c=1:nCores
		coreData.Nodes.Ghosts_ImReceiving{c} = sort(unique(coreData.Nodes.Ghosts_ImReceiving{c}));
		coreData.Nodes.Ghosts_ImSending{c} = sort(unique(coreData.Nodes.Ghosts_ImSending{c}));
	end

	% renumbering
	[~,reorderNodes] = sort(coreData.NodeLocs);
	Nodes = Nodes(reorderNodes,:);
	coreData.NodeLocs = coreData.NodeLocs(reorderNodes);

	for ng=1:length(NodeGroup)
		[~,reorderNodeGroup{ng}] = sort(coreData.NodeGroup{ng}.core);
		coreData.NodeGroup{ng}.core = coreData.NodeGroup{ng}.core(reorderNodeGroup{ng});
		NOLD = NodeGroup{ng}.Nodes;
		for n=1:length(NOLD)
			ToRenumber = NOLD(reorderNodeGroup{ng}(n));
			Renumbered = find(ToRenumber==reorderNodes);
			NodeGroup{ng}.Nodes(n) = Renumbered;
		end
	end

	for eg=1:length(ElementGroup)
		[~,reorderElems{eg}] = sort(coreData.ElementGroup{eg}.core);
		coreData.ElementGroup{eg}.core = coreData.ElementGroup{eg}.core(reorderElems{eg});
		EOLD = ElementGroup{eg}.Elements;

		for n=1:size(EOLD,1)
			for nn=1:size(EOLD,2)
				ToRenumber = EOLD(reorderElems{eg}(n),nn);
				Renumbered = find(ToRenumber==reorderNodes);
				ElementGroup{eg}.Elements(n,nn) = Renumbered;
			end
		end
	end

	for c=1:nCores
		for n=1:length(coreData.Nodes.PrivateNodes{c})
			coreData.Nodes.PrivateNodes{c}(n) = find(coreData.Nodes.PrivateNodes{c}(n)==reorderNodes);
		end
		for n=1:length(coreData.Nodes.Ghosts_ImReceiving{c})
			coreData.Nodes.Ghosts_ImReceiving{c}(n) = find(coreData.Nodes.Ghosts_ImReceiving{c}(n)==reorderNodes);
		end
		for n=1:length(coreData.Nodes.Ghosts_ImSending{c})
			coreData.Nodes.Ghosts_ImSending{c}(n) = find(coreData.Nodes.Ghosts_ImSending{c}(n)==reorderNodes);
		end
		coreData.Nodes.PrivateNodes{c} = sort(coreData.Nodes.PrivateNodes{c});
		coreData.Nodes.Ghosts_ImReceiving{c} = sort(coreData.Nodes.Ghosts_ImReceiving{c});
		coreData.Nodes.Ghosts_ImSending{c} = sort(coreData.Nodes.Ghosts_ImSending{c});
	end
	
	for c=1:nCores
		NodeRange = find(coreData.NodeLocs==c);
		coreData.ToSave{c}.Noderange = uint64([min(NodeRange) max(NodeRange)]);
		coreData.ToSave{c}.Ghosts = uint64(coreData.Nodes.Ghosts_ImReceiving{c});

		for ng=1:length(NodeGroup)
			myNodes = find(coreData.NodeGroup{ng}.core==c);
			if isempty(myNodes)
				coreData.ToSave{c}.hasnodegroup(ng) = false;
				coreData.ToSave{c}.NodeGrouprange(ng,:) = [nan nan];
			else
				coreData.ToSave{c}.hasnodegroup(ng) = true;
				coreData.ToSave{c}.NodeGrouprange(ng,:) = uint64([min(myNodes) max(myNodes)]);
			end
		end

		for eg=1:length(coreData.ElementGroup)
			myElems = find(coreData.ElementGroup{eg}.core==c);
			if isempty(myElems)
				coreData.ToSave{c}.haselems(eg) = false;
				coreData.ToSave{c}.Elemrange(eg,:) = [nan nan];
			else
				coreData.ToSave{c}.haselems(eg) = true;
				coreData.ToSave{c}.Elemrange(eg,:) = uint64([min(myElems) max(myElems)]);
			end
		end

	end
end





function Elementgroups = getBoundaryGroup(Elementgroups, GroupIndex, model, Face)
	Surfaces = [1 2 3 5 6 7;
		        1 2 4 5 9 8;
				1 3 4 7 10 8;
				2 3 4 6 10 9];

	Elementgroups{GroupIndex}.Elements = [];
	
	BoundaryNodes = findNodes(model,'region','Face',Face);
	iElem = 0;
	for i=1:length(Elementgroups{1}.Elements)
		V_Elem = Elementgroups{1}.Elements(i,:);
		
		matchedNodes = [];
		for j=1:length(V_Elem)
			ix = find(V_Elem(j)==BoundaryNodes);
			if (length(ix)>0)
				matchedNodes(end+1) = j;

			end
		end
		if (length(matchedNodes) == 6) %element on interface STILL REQUIRES SORTING
			SurfFound = -1;
			for n=1:4
				tf = isequal(sort(matchedNodes), sort(Surfaces(n,:)));
				if (tf)
					SurfFound = n;
				end
			end
			Elementgroups{GroupIndex}.Elements(end+1,:) = V_Elem(Surfaces(SurfFound,:));
		end
	end
end

function PlotMesh(Nodes,Elementgroups,groupsToPlot)
	cArray = ["r", "g", "b", "c", "y", "m"];
	for g=groupsToPlot
		F = Elementgroups{g}.Elements(:,[1 2 3]);
		patch('Faces',F,'Vertices',Nodes,'FaceColor',cArray(g),'EdgeColor','k','FaceAlpha',.5)
	end
end


function saveToHDF(fname, Nodes, NodeGroup, ElementGroup, coreData)
	if exist(fname+".h5", 'file')==2
		delete(fname+".h5");
	end

	% nodes
	h5create(fname+".h5",'/nodes',size(Nodes),'Datatype','double');
	h5write(fname+".h5",'/nodes', Nodes)

	% nodegroups
	h5create(fname+".h5",'/nodegroupnames',length(NodeGroup),'Datatype','string');
	for i=1:length(NodeGroup)
		names{i} = NodeGroup{i}.Name;
	end

	h5write(fname+".h5",'/nodegroupnames', string(names))
	for i=1:length(NodeGroup)
		h5create(fname+".h5",'/nodegroups/'+NodeGroup{i}.Name,length(NodeGroup{i}.Nodes),'Datatype','uint64');
		h5write(fname+".h5",'/nodegroups/'+NodeGroup{i}.Name,NodeGroup{i}.Nodes-1);
	end

	% Elementgroups
	clear names types
	h5create(fname+".h5",'/elementgroupnames',length(ElementGroup),'Datatype','string');
	h5create(fname+".h5",'/elementgrouptypes',length(ElementGroup),'Datatype','string');
	for i=1:length(ElementGroup)
		names{i} = ElementGroup{i}.Name;
		types{i} = ElementGroup{i}.Type;
	end

	h5write(fname+".h5",'/elementgroupnames', string(names))
	h5write(fname+".h5",'/elementgrouptypes', string(types))
	for i=1:length(ElementGroup)
		h5create(fname+".h5",'/elementgroups/'+ElementGroup{i}.Name,size(ElementGroup{i}.Elements),'Datatype','uint64');
		h5write(fname+".h5",'/elementgroups/'+ElementGroup{i}.Name,ElementGroup{i}.Elements-1);
	end

	% partition data
	h5create(fname+".h5",'/NCores', 1,'Datatype','uint64');
	h5write(fname+".h5",'/NCores', uint64(coreData.NCores));

	for c=1:coreData.NCores
		h5create(fname+".h5",'/Partition/'+string(c-1)+'/noderange', 2,'Datatype','uint64');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/noderange',coreData.ToSave{c}.Noderange-1);

		if (isempty(coreData.ToSave{c}.Ghosts))
			h5create(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', 1,'Datatype','uint8');
			h5write(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', uint8(false));
		else
			h5create(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', 1,'Datatype','uint8');
			h5write(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', uint8(true));

			h5create(fname+".h5",'/Partition/'+string(c-1)+'/ghosts', length(coreData.ToSave{c}.Ghosts),'Datatype','uint64');
			h5write(fname+".h5",'/Partition/'+string(c-1)+'/ghosts',coreData.ToSave{c}.Ghosts-1);
		end


		h5create(fname+".h5",'/Partition/'+string(c-1)+'/Haselems', length(coreData.ToSave{c}.haselems),'Datatype','uint8');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/Haselems',uint8(coreData.ToSave{c}.haselems));

		h5create(fname+".h5",'/Partition/'+string(c-1)+'/elemrange', size(coreData.ToSave{c}.Elemrange),'Datatype','uint64');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/elemrange',coreData.ToSave{c}.Elemrange-1);

		h5create(fname+".h5",'/Partition/'+string(c-1)+'/hasnodegroup', length(coreData.ToSave{c}.hasnodegroup),'Datatype','uint8');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/hasnodegroup',uint8(coreData.ToSave{c}.hasnodegroup));

		h5create(fname+".h5",'/Partition/'+string(c-1)+'/nodegrouprange', size(coreData.ToSave{c}.NodeGrouprange),'Datatype','uint64');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/nodegrouprange',coreData.ToSave{c}.NodeGrouprange-1);
	end

	%h5disp(fname+".h5")
end
