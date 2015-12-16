
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
function showChain(){
    var i=prompt("",window.show_val)
    if(!i&&i!==0) return;
    window.show_val=i
	$.each(a.coords,function(){
		a.markers[this.index].marker.setIcon(OSRM.G.icons['marker-drag']);
		a.markers[this.index].marker.setLabel(this.index)
	})
	if(i!="00")
		$.each(JSON.parse("["+i+"]"),function(){
			a.markers[this].marker.setIcon(OSRM.G.icons['marker-via']);
			a.markers[this].marker.setLabel(this)
		})
}

function lineChain(from, to)
{ return [a.markers[from].position, a.markers[to].position] }



var tails=new OSRM.MultiRoute("tails");
var inners=new OSRM.MultiRoute("inners");
var a={coords:[]}

function SCT(b){
$.each(a.coords,function(){
	a.markers[this.index].hide()
})
a=b;
a.markers=[];
$.each(a.coords,function(){
	var tit="";
	if(this.tws)
    $.each(this.tws,function(){
		var start = new Date(this.start*1000);
		var finish = new Date(this.finish*1000);
		var sh = start.getUTCHours();
		var sm = "0" + start.getUTCMinutes();
		var fh = finish.getUTCHours();
		var fm = "0" + finish.getUTCMinutes();

		tit+=sh+ ':' + sm.substr(sm.length-2)+ ' - ' + fh+ ':' + fm.substr(fm.length-2)+"\n"
	})
	a.markers[this.index]=new OSRM.RouteMarker(OSRM.C.VIA_LABEL, { draggable: 0, icon: OSRM.G.icons["marker-via"], dragicon: OSRM.G.icons["marker-via-drag"], title:tit}, this)
	a.markers[this.index].show();
	a.markers[this.index].marker.setIcon(OSRM.G.icons['marker-drag']);
	a.markers[this.index].marker.setLabel(this.index)
})

$.each(a.cores,function(){
	var core = this;
	$.each(core.ptr_wrapper.data.inners,function(){
		a.markers[this].host_core=core.ptr_wrapper.data.index;
	})
})

tails.clearRoutes()
inners.clearRoutes()
$.each(a.nearest_graph.forward_data,function(){
	var from=this.key;
	$.each(this.value,function(){
		if(a.markers[from].host_core==a.markers[this].host_core && a.markers[from].host_core!=undefined)
			inners.addRoute(lineChain(from, this));
		else
			tails.addRoute(lineChain(from, this));
	})
})
tails.show();
inners.show();
inners.setStyle({color:'#222222', weight:2, dashArray:"8,6"});
tails.setStyle({color:'#3366FF', weight:1, dashArray:""});

var box=document.getElementById('clusters-box');
var text='<a href="#" onclick="return showAll()">showAll</a> <a href="#" onclick="return showChain()">showChain</a><br/><br/>';
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
}