a.markers=[];$.each(a.coords,function(){
	a.markers[this.index]=new OSRM.RouteMarker(OSRM.C.VIA_LABEL, { draggable: !0, icon: OSRM.G.icons["marker-via"], dragicon: OSRM.G.icons["marker-via-drag"]}, this)
	a.markers[this.index].show();
})
$.each(a.coords,function(){
	a.markers[this.index].marker.setIcon(OSRM.G.icons['marker-drag']);
	a.markers[this.index].marker.setLabel(this.index)
})
function showAll(){
	$.each(a.coords,function(){
		a.markers[this.index].marker.setOpacity(1);
	})
}
function showCore(c){
	$.each(a.coords,function(){
		a.markers[this.index].marker.setOpacity(0.1);
	})
	$.each(a.cores[c].ptr_wrapper.data.tails,function(){
		a.markers[this].marker.setOpacity(0.5);
	})
	$.each(a.cores[c].ptr_wrapper.data.inners,function(){
		a.markers[this].marker.setOpacity(1);
	})
}

function lineChain(from, to)
{ return [a.markers[from].position, a.markers[to].position] }



$.each(a.cores,function(){
	var core = this;
	$.each(core.ptr_wrapper.data.inners,function(){
		a.markers[this].host_core=core.ptr_wrapper.data.index;
	})
})

var tails=new OSRM.MultiRoute("tails");
var inners=new OSRM.MultiRoute("inners");
$.each(a.nearest_graph.forward_data,function(){
	var from=this.key;
	$.each(this.value,function(){
		if(a.markers[from].host_core==a.markers[this].host_core && a.markers[from].host_core!=undefined)
			tails.addRoute(lineChain(from, this));
		else
			inners.addRoute(lineChain(from, this));
	})
})
tails.show();
inners.show();
tails.setStyle({color:'#222222', weight:2, dashArray:"8,6"});
inners.setStyle({color:'#3366FF', weight:1, dashArray:""});

var box=document.getElementById('clusters-box');
var text='<a href="#" onclick="return showAll()">showAll</a><br/><br/>';
text+='Ядра<ul>';
$.each(a.cores,function(i){
	var core = this;
	text+='<li><a href="#" onclick="return showCore('+i+')">'+core.ptr_wrapper.data.index+'</a>='+core.ptr_wrapper.data.inners[0]
	for(var j = 1; j<core.ptr_wrapper.data.inners.length; ++j)
		text+='+'+core.ptr_wrapper.data.inners[j]
	text+='</li>';
})
text+='</ul>';
box.innerHTML=text;