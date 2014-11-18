L.DomUtil.changeClass = function(a, b, c) {
    a.className = L.Util.trim((" " + a.className + " ").replace(" " + b + " ", " " + c + " "))
}, L.LabelMarker = L.Marker.extend({
    changeIcon: function(a) {
        this.options.icon = a, this._map && this._changeIcon()
    },
    setLabel: function(a) {
        this._icon && (this._icon.lastChild.innerHTML = a, this._icon.lastChild.style.display = "block")
    },
    setTitle: function(a) {
        this.options.title = a, this._icon.title = a
    },
    _changeIcon: function() {
        var a = this.options,
            b = this._map,
            c = b.options.zoomAnimation && b.options.markerZoomAnimation,
            d = c ? "leaflet-zoom-animated" : "leaflet-zoom-hide";
        this._icon && (this._icon = a.icon.changeIcon(this._icon), L.DomUtil.addClass(this._icon, d), L.DomUtil.addClass(this._icon, "leaflet-clickable"))
    }
}), L.LabelMarkerIcon = L.Icon.extend({
    _createImg: function(a) {
        var b;
        if (L.Browser.ie6) b = document.createElement("div"), b.style.filter = 'progid:DXImageTransform.Microsoft.AlphaImageLoader(src="' + a + '")';
        else {
            b = document.createElement("div");
            var c = document.createElement("img"),
                d = document.createElement("div");
            c.src = a, d.className = "via-counter", d.innerHTML = "", b.appendChild(c), b.appendChild(d)
        }
        return b
    },
    changeIcon: function(a) {
        return this._changeIcon("icon", a)
    },
    changeShadow: function(a) {
        return this.options.shadowUrl ? this._changeIcon("shadow", a) : null
    },
    _changeIcon: function(a, b) {
        var c = this._getIconUrl(a);
        if (!c) {
            if ("icon" === a) throw new Error("iconUrl not set in Icon options (see the docs).");
            return null
        }
        var d = this._changeImg(c, b);
        return this._setIconStyles(d, a), d
    },
    _changeImg: function(a, b) {
        return L.Browser.ie6 ? b.style.filter = 'progid:DXImageTransform.Microsoft.AlphaImageLoader(src="' + a + '")' : b.firstChild.src = a, b
    }
});
var OSRM = {};
OSRM.VERSION = "0.1.11", OSRM.DATE = "131122", OSRM.CONSTANTS = {}, OSRM.DEFAULTS = {}, OSRM.GLOBALS = {}, OSRM.Control = {}, OSRM.G = OSRM.GLOBALS, OSRM.C = OSRM.CONSTANTS, OSRM.DEFAULTS = {
        ROUTING_ENGINES: [{
            url: "http://192.168.9.29:5000/viaroute",
            timestamp: "http://192.168.9.29:5000/timestamp",
            table: "http://192.168.9.29:5000/table",
            metric: 0,
            label: "ENGINE_0"
        }],
        WEBSITE_URL: document.URL.replace(/#*(\?.*|$)/i, ""),
        HOST_GEOCODER_URL: "https://nominatim.openstreetmap.org/search",
        HOST_REVERSE_GEOCODER_URL: "https://nominatim.openstreetmap.org/reverse",
        HOST_SHORTENER_URL: "http://short.project-osrm.org/",
        SHORTENER_PARAMETERS: "%url&jsonp=%jsonp",
        SHORTENER_REPLY_PARAMETER: "ShortURL",
        ROUTING_ENGINE: 0,
        DISTANCE_FORMAT: 0,
        GEOCODER_BOUNDS: "",
        ZOOM_LEVEL: 14,
        HIGHLIGHT_ZOOM_LEVEL: 16,
        JSONP_TIMEOUT: 1e4,
        EDITOR_MIN_ZOOM_LEVEL: 16,
        JOSM_MIN_ZOOM_LEVEL: 16,
        NOTES_MIN_ZOOM_LEVEL: 8,
        ONLOAD_ZOOM_LEVEL: 5,
        ONLOAD_LATITUDE: 48.84,
        ONLOAD_LONGITUDE: 10.1,
        ONLOAD_SOURCE: "",
        ONLOAD_TARGET: "",
        LANGUAGE: "en",
        LANGUAGE_USE_BROWSER_SETTING: !0,
        LANUGAGE_ONDEMAND_RELOADING: !0,
        LANGUAGE_SUPPORTED: [{
            encoding: "en",
            name: "English"
        }, {
            encoding: "bg",
            name: "Български"
        }, {
            encoding: "ca",
            name: "Català"
        }, {
            encoding: "cs",
            name: "Česky"
        }, {
            encoding: "de",
            name: "Deutsch"
        }, {
            encoding: "da",
            name: "Dansk"
        }, {
            encoding: "el",
            name: "Ελληνικά"
        }, {
            encoding: "eo",
            name: "Esperanto"
        }, {
            encoding: "es",
            name: "Español"
        }, {
            encoding: "fi",
            name: "Suomi"
        }, {
            encoding: "fr",
            name: "Français"
        }, {
            encoding: "it",
            name: "Italiano"
        }, {
            encoding: "ja",
            name: "日本語"
        }, {
            encoding: "ka",
            name: "ქართული"
        }, {
            encoding: "lv",
            name: "Latviešu"
        }, {
            encoding: "nb",
            name: "Bokmål"
        }, {
            encoding: "nl",
            name: "Nederlands"
        }, {
            encoding: "pl",
            name: "Polski"
        }, {
            encoding: "pt",
            name: "Português"
        }, {
            encoding: "ro",
            name: "Română"
        }, {
            encoding: "ru",
            name: "Русский"
        }, {
            encoding: "sk",
            name: "Slovensky"
        }, {
            encoding: "sv",
            name: "Svenska"
        }, {
            encoding: "ta",
            name: "தமிழ்"
        }, {
            encoding: "tr",
            name: "Türkçe"
        }, {
            encoding: "uk",
            name: "Українська"
        }],
        TILE_SERVERS: [{
            display_name: "OSRM bright",
            url: "https://{s}.tiles.mapbox.com/v4/dennisl.4e2aab76/{z}/{x}/{y}.png?access_token=pk.eyJ1IjoiZGVubmlzbCIsImEiOiI4bkR5ZDNJIn0.KuIIMvojByGZzIh1RCzKwg",
            attribution: '<a href="https://www.mapbox.com/about/maps">© Mapbox</a> <a href="http://openstreetmap.org/copyright">© OpenStreetMap</a> | <a href="http://mapbox.com/map-feedback/">Improve this map</a>',
            options: {
                maxZoom: 18
            }
        }, {
            display_name: "Mapbox Terrain",
            url: "https://{s}.tiles.mapbox.com/v3/dennisl.map-dfbkqsr2/{z}/{x}/{y}.png",
            attribution: '<a href="https://www.mapbox.com/about/maps">© Mapbox</a> <a href="http://openstreetmap.org/copyright">© OpenStreetMap</a> | <a href="http://mapbox.com/map-feedback/">Improve this map</a>',
            options: {
                maxZoom: 18
            }
        }, {
            display_name: "Mapbox Labelled Satellite",
            url: "https://{s}.tiles.mapbox.com/v3/dennisl.map-6g3jtnzm/{z}/{x}/{y}.png",
            attribution: '<a href="https://www.mapbox.com/about/maps">© Mapbox</a> <a href="http://openstreetmap.org/copyright">© OpenStreetMap</a> | <a href="http://mapbox.com/map-feedback/">Improve this map</a>',
            options: {
                maxZoom: 18
            }
        }, {
            display_name: "Mapbox Satellite",
            url: "https://{s}.tiles.mapbox.com/v3/dennisl.map-inp5al1s/{z}/{x}/{y}.png",
            attribution: '<a href="https://www.mapbox.com/about/maps">© Mapbox</a> <a href="http://openstreetmap.org/copyright">© OpenStreetMap</a> | <a href="http://mapbox.com/map-feedback/">Improve this map</a>',
            options: {
                maxZoom: 18
            }
        }, {
            display_name: "osm.org",
            url: "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
            attribution: '© <a href="http://www.openstreetmap.org/copyright/en">OpenStreetMap</a> contributors',
            options: {
                maxZoom: 18
            }
        }, {
            display_name: "osm.de",
            url: "http://{s}.tile.openstreetmap.de/tiles/osmde/{z}/{x}/{y}.png",
            attribution: '© <a href="http://www.openstreetmap.org/copyright/en">OpenStreetMap</a> contributors',
            options: {
                maxZoom: 18
            }
        }, {
            display_name: "MapQuest",
            url: "http://otile{s}.mqcdn.com/tiles/1.0.0/osm/{z}/{x}/{y}.png",
            attribution: '© <a href="http://www.openstreetmap.org/copyright/en">OpenStreetMap</a> contributors, Imagery © <a href="http://www.mapquest.de/">MapQuest</a>',
            options: {
                maxZoom: 18,
                subdomains: "1234"
            }
        }],
        OVERLAY_SERVERS: [{
            display_name: "Small Components",
            url: "http://tools.geofabrik.de/osmi/tiles/routing_i/{z}/{x}/{y}.png",
            attribution: "",
            options: {}
        }],
        NOTIFICATIONS: {
            LOCALIZATION: 18e5,
            CLICKING: 6e4,
            DRAGGING: 12e4,
            MAINTENANCE: !1
        },
        OVERRIDE_MAINTENANCE_NOTIFICATION_HEADER: void 0,
        OVERRIDE_MAINTENANCE_NOTIFICATION_BODY: void 0
    },
    function() {
        var a = navigator.userAgent;
        OSRM.Browser = {
            FF3: a.search(/Firefox\/3/),
            IE6_7: a.search(/MSIE (6|7)/),
            IE6_8: a.search(/MSIE (6|7|8)/),
            IE6_9: a.search(/MSIE (6|7|8|9)/)
        }
    }(), document.head = document.head || document.getElementsByTagName("head")[0], OSRM.Browser.getElementsByClassName = function(a, b) {
        for (var c = [], d = new RegExp("(^| )" + b + "( |$)"), e = a.getElementsByTagName("*"), f = 0, g = e.length; g > f; f++) d.test(e[f].className) && c.push(e[f]);
        return c
    }, OSRM.Browser.onLoadHandler = function(a, b) {
        b = b || window;
        var c = b.document;
        if (b.addEventListener) {
            var d = function() {
                b.removeEventListener("DOMContentLoaded", arguments.callee, !1), a.call()
            };
            b.addEventListener("DOMContentLoaded", d, !1)
        } else if (c.attachEvent) {
            var d = function() {
                ("interactive" === c.readyState || "complete" === c.readyState) && (c.detachEvent("onreadystatechange", arguments.callee), a.call())
            };
            c.attachEvent("onreadystatechange", d)
        }
    }, OSRM.Browser.onUnloadHandler = function(a, b) {
        b = b || window;
        var c = b.document;
        b.addEventListener ? b.addEventListener("unload", a, !1) : c.attachEvent && c.attachEvent("onunload", a)
    },
    function() {
        var a = function() {};
        OSRM.inheritFrom = function(b, c) {
            a.prototype = c.prototype, b.prototype = new a, b.prototype.constructor = b, b.prototype.base = c.prototype
        }
    }(), OSRM.extend = function(a, b) {
        for (property in b) a.prototype[property] = b[property]
    }, OSRM.bind = function(a, b) {
        return function() {
            b.apply(a, arguments)
        }
    }, OSRM.concat = function(a, b) {
        return function() {
            a.apply(this, arguments), b.apply(this, arguments)
        }
    }, OSRM.init = function() {
        1 != OSRM.checkOldBrowser() && (OSRM.showHTML(), OSRM.prefetchImages(), OSRM.prefetchIcons(), OSRM.prefetchCSSIcons(), OSRM.GUI.init(), OSRM.Map.init(), OSRM.Printing.init(), OSRM.Routing.init(), OSRM.RoutingAlternatives.init(), OSRM.Localization.init(), OSRM.parseParameters(), 1 != OSRM.GUI.inMaintenance() && (0 == OSRM.G.initial_position_override && OSRM.Map.initPosition(), OSRM.Map.initFinally()))
    }, OSRM.GLOBALS.images = {}, OSRM.prefetchImages = function() {
        for (var a = [{
                id: "marker-shadow",
                url: "leaflet/images/marker-shadow.png"
            }, {
                id: "marker-source",
                url: "images/marker-source.png"
            }, {
                id: "marker-target",
                url: "images/marker-target.png"
            }, {
                id: "marker-via",
                url: "images/marker-via.png"
            }, {
                id: "marker-highlight",
                url: "images/marker-highlight.png"
            }, {
                id: "marker-source-drag",
                url: "images/marker-source-drag.png"
            }, {
                id: "marker-target-drag",
                url: "images/marker-target-drag.png"
            }, {
                id: "marker-via-drag",
                url: "images/marker-via-drag.png"
            }, {
                id: "marker-highlight-drag",
                url: "images/marker-highlight-drag.png"
            }, {
                id: "marker-drag",
                url: "images/marker-drag.png"
            }, {
                id: "cancel",
                url: "images/cancel.png"
            }, {
                id: "cancel_active",
                url: "images/cancel_active.png"
            }, {
                id: "cancel_hover",
                url: "images/cancel_hover.png"
            }, {
                id: "restore",
                url: "images/restore.png"
            }, {
                id: "restore_active",
                url: "images/restore_active.png"
            }, {
                id: "restore_hover",
                url: "images/restore_hover.png"
            }, {
                id: "up",
                url: "images/up.png"
            }, {
                id: "up_active",
                url: "images/up_active.png"
            }, {
                id: "up_hover",
                url: "images/up_hover.png"
            }, {
                id: "down",
                url: "images/down.png"
            }, {
                id: "down_active",
                url: "images/down_active.png"
            }, {
                id: "down_hover",
                url: "images/down_hover.png"
            }, {
                id: "config",
                url: "images/config.png"
            }, {
                id: "config_active",
                url: "images/config_active.png"
            }, {
                id: "config_hover",
                url: "images/config_hover.png"
            }, {
                id: "mapping",
                url: "images/mapping.png"
            }, {
                id: "mapping_active",
                url: "images/mapping_active.png"
            }, {
                id: "mapping_hover",
                url: "images/mapping_hover.png"
            }, {
                id: "printer",
                url: "images/printer.png"
            }, {
                id: "printer_active",
                url: "images/printer_active.png"
            }, {
                id: "printer_hover",
                url: "images/printer_hover.png"
            }, {
                id: "printer_inactive",
                url: "images/printer_inactive.png"
            }, {
                id: "zoom_in",
                url: "images/zoom_in.png"
            }, {
                id: "zoom_in_active",
                url: "images/zoom_in_active.png"
            }, {
                id: "zoom_in_hover",
                url: "images/zoom_in_hover.png"
            }, {
                id: "zoom_in_inactive",
                url: "images/zoom_in_inactive.png"
            }, {
                id: "zoom_out",
                url: "images/zoom_out.png"
            }, {
                id: "zoom_out_active",
                url: "images/zoom_out_active.png"
            }, {
                id: "zoom_out_hover",
                url: "images/zoom_out_hover.png"
            }, {
                id: "zoom_out_inactive",
                url: "images/zoom_out_inactive.png"
            }, {
                id: "locations_user",
                url: "images/locations_user.png"
            }, {
                id: "locations_user_active",
                url: "images/locations_user_active.png"
            }, {
                id: "locations_user_hover",
                url: "images/locations_user_hover.png"
            }, {
                id: "locations_user_inactive",
                url: "images/locations_user_inactive.png"
            }, {
                id: "locations_route",
                url: "images/locations_route.png"
            }, {
                id: "locations_route_active",
                url: "images/locations_route_active.png"
            }, {
                id: "locations_route_hover",
                url: "images/locations_route_hover.png"
            }, {
                id: "locations_route_inactive",
                url: "images/locations_route_inactive.png"
            }, {
                id: "layers",
                url: "images/layers.png"
            }, {
                id: "direction_0",
                url: "images/default.png"
            }, {
                id: "direction_1",
                url: "images/continue.png"
            }, {
                id: "direction_2",
                url: "images/slight-right.png"
            }, {
                id: "direction_3",
                url: "images/turn-right.png"
            }, {
                id: "direction_4",
                url: "images/sharp-right.png"
            }, {
                id: "direction_5",
                url: "images/u-turn.png"
            }, {
                id: "direction_6",
                url: "images/sharp-left.png"
            }, {
                id: "direction_7",
                url: "images/turn-left.png"
            }, {
                id: "direction_8",
                url: "images/slight-left.png"
            }, {
                id: "direction_10",
                url: "images/head.png"
            }, {
                id: "direction_11",
                url: "images/round-about.png"
            }, {
                id: "direction_15",
                url: "images/target.png"
            }, {
                id: "osrm-logo",
                url: "images/osrm-logo.png"
            }, {
                id: "selector",
                url: "images/selector.png"
            }], b = 0; b < a.length; b++) OSRM.G.images[a[b].id] = new Image, OSRM.G.images[a[b].id].src = a[b].url
    }, OSRM.GLOBALS.icons = {}, OSRM.prefetchIcons = function() {
        for (var a = [{
                id: "marker-source",
                image_id: "marker-source"
            }, {
                id: "marker-target",
                image_id: "marker-target"
            }, {
                id: "marker-via",
                image_id: "marker-via"
            }, {
                id: "marker-highlight",
                image_id: "marker-highlight"
            }, {
                id: "marker-source-drag",
                image_id: "marker-source-drag"
            }, {
                id: "marker-target-drag",
                image_id: "marker-target-drag"
            }, {
                id: "marker-via-drag",
                image_id: "marker-via-drag"
            }, {
                id: "marker-highlight-drag",
                image_id: "marker-highlight-drag"
            }], b = L.LabelMarkerIcon.extend({
                options: {
                    shadowUrl: OSRM.G.images["marker-shadow"].getAttribute("src"),
                    iconSize: [25, 41],
                    shadowSize: [41, 41],
                    iconAnchor: [13, 41],
                    shadowAnchor: [13, 41],
                    popupAnchor: [0, -33]
                }
            }), c = 0; c < a.length; c++) OSRM.G.icons[a[c].id] = new b({
            iconUrl: OSRM.G.images[a[c].image_id].getAttribute("src")
        });
        OSRM.G.icons["marker-drag"] = new L.LabelMarkerIcon({
            iconUrl: OSRM.G.images["marker-drag"].getAttribute("src"),
            iconSize: new L.Point(18, 18)
        })
    }, OSRM.prefetchCSSIcons = function() {
        for (var a = [{
                id: "#gui-printer-inactive",
                image_id: "printer_inactive"
            }, {
                id: "#gui-printer",
                image_id: "printer"
            }, {
                id: "#gui-printer:hover",
                image_id: "printer_hover"
            }, {
                id: "#gui-printer:active",
                image_id: "printer_active"
            }, {
                id: ".gui-zoom-in-inactive",
                image_id: "zoom_in_inactive"
            }, {
                id: ".gui-zoom-in",
                image_id: "zoom_in"
            }, {
                id: ".gui-zoom-in:hover",
                image_id: "zoom_in_hover"
            }, {
                id: ".gui-zoom-in:active",
                image_id: "zoom_in_active"
            }, {
                id: ".gui-zoom-out-inactive",
                image_id: "zoom_out_inactive"
            }, {
                id: ".gui-zoom-out",
                image_id: "zoom_out"
            }, {
                id: ".gui-zoom-out:hover",
                image_id: "zoom_out_hover"
            }, {
                id: ".gui-zoom-out:active",
                image_id: "zoom_out_active"
            }, {
                id: ".gui-locations-user-inactive",
                image_id: "locations_user_inactive"
            }, {
                id: ".gui-locations-user",
                image_id: "locations_user"
            }, {
                id: ".gui-locations-user:hover",
                image_id: "locations_user_hover"
            }, {
                id: ".gui-locations-user:active",
                image_id: "locations_user_active"
            }, {
                id: ".gui-locations-route-inactive",
                image_id: "locations_route_inactive"
            }, {
                id: ".gui-locations-route",
                image_id: "locations_route"
            }, {
                id: ".gui-locations-route:hover",
                image_id: "locations_route_hover"
            }, {
                id: ".gui-locations-route:active",
                image_id: "locations_route_active"
            }, {
                id: ".gui-layers",
                image_id: "layers"
            }, {
                id: ".cancel-marker",
                image_id: "cancel"
            }, {
                id: ".cancel-marker:hover",
                image_id: "cancel_hover"
            }, {
                id: ".cancel-marker:active",
                image_id: "cancel_active"
            }, {
                id: ".up-marker",
                image_id: "up"
            }, {
                id: ".up-marker:hover",
                image_id: "up_hover"
            }, {
                id: ".up-marker:active",
                image_id: "up_active"
            }, {
                id: ".down-marker",
                image_id: "down"
            }, {
                id: ".down-marker:hover",
                image_id: "down_hover"
            }, {
                id: ".down-marker:active",
                image_id: "down_active"
            }, {
                id: "#input-mask-header",
                image_id: "osrm-logo"
            }, {
                id: ".styled-select",
                image_id: "selector"
            }, {
                id: "#config-handle-icon",
                image_id: "config"
            }, {
                id: "#config-handle-icon:hover",
                image_id: "config_hover"
            }, {
                id: "#config-handle-icon:active",
                image_id: "config_active"
            }, {
                id: "#mapping-handle-icon",
                image_id: "mapping"
            }, {
                id: "#mapping-handle-icon:hover",
                image_id: "mapping_hover"
            }, {
                id: "#mapping-handle-icon:active",
                image_id: "mapping_active"
            }, {
                id: "#main-handle-icon",
                image_id: "restore"
            }, {
                id: "#main-handle-icon:hover",
                image_id: "restore_hover"
            }, {
                id: "#main-handle-icon:active",
                image_id: "restore_active"
            }], b = OSRM.CSS.getStylesheet("main.css"), c = 0; c < a.length; c++) OSRM.CSS.insert(b, a[c].id, "background-image:url(" + OSRM.G.images[a[c].image_id].getAttribute("src") + ");")
    }, OSRM.checkOldBrowser = function() {
        return -1 == OSRM.Browser.IE6_7 ? !1 : (document.getElementById("old-browser-warning").style.display = "block", !0)
    }, OSRM.showHTML = function() {
        document.getElementById("map").style.display = "block", document.getElementById("gui").style.display = "block"
    }, OSRM.parseParameters = function() {
        var a = document.location.search.substr(1, document.location.search.length);
        if (OSRM.G.initial_position_override = !1, !(a.length > 1e3 || 0 == a.length)) {
            var b = {};
            b.active_routing_engine = OSRM.DEFAULTS.ROUTING_ENGINE;
            for (var c = a.split("&"), d = 0; d < c.length; d++) {
                var e = c[d].split("=");
                if (2 == e.length)
                    if ("hl" == e[0]) OSRM.Localization.setLanguage(e[1]);
                    else if ("df" == e[0]) {
                    var f = parseInt(e[1]);
                    if (0 != f && 1 != f) return;
                    OSRM.GUI.setDistanceFormat(f)
                } else if ("loc" == e[0]) b.locations = b.locations || [], b.locations.push(decodeURI(e[1]).replace(/<\/?[^>]+(>|$)/g, ""));
                else if ("dest" == e[0]) b.destinations = b.dlocations || [], b.destinations.push(decodeURI(e[1]).replace(/<\/?[^>]+(>|$)/g, ""));
                else if ("destname" == e[0]) b.destination_name = decodeURI(e[1]).replace(/<\/?[^>]+(>|$)/g, "");
                else if ("z" == e[0]) {
                    var g = Number(e[1]);
                    if (0 > g || g > 18) return;
                    b.zoom = g
                } else if ("center" == e[0]) {
                    var h = unescape(e[1]).split(",");
                    if (2 != h.length || !OSRM.Utils.isLatitude(h[0]) || !OSRM.Utils.isLongitude(h[1])) return;
                    b.center = new L.LatLng(h[0], h[1])
                } else if ("alt" == e[0]) {
                    var i = Number(e[1]);
                    if (0 > i || i > OSRM.RoutingAlternatives > 10) return;
                    b.active_alternative = i
                } else if ("re" == e[0]) {
                    var j = Number(e[1]);
                    if (0 > j || j >= OSRM.DEFAULTS.ROUTING_ENGINES.length) return;
                    b.active_routing_engine = j
                } else if ("ly" == e[0]) {
                    for (var k = Number(e[1]), l = 0; l < OSRM.DEFAULTS.TILE_SERVERS.length; l++)
                        if (OSRM.Utils.getHash(OSRM.DEFAULTS.TILE_SERVERS[l].display_name) == k) {
                            OSRM.G.map.layerControl.setActiveLayerByName(OSRM.DEFAULTS.TILE_SERVERS[l].display_name);
                            break
                        }
                } else "notifications" == e[0] ? "hide" == e[1] && OSRM.GUI.deactivateTooltips() : "mainbox" == e[0] && "hide" == e[1] && (OSRM.G.main_handle.$hideBox(), OSRM.G.map.zoomControl.show())
            }
            if (1 != OSRM.GUI.inMaintenance())
                if (OSRM.GUI.setRoutingEngine(b.active_routing_engine), b.locations || b.destinations) {
                    var m = b.destinations ? b.destinations : b.locations,
                        n = b.destinations ? "_showInitResults_Destinations" : "_showInitResults_Locations";
                    OSRM.G.initial_positions = {};
                    var o = OSRM.G.initial_positions;
                    o.positions = [], o.done = 0, o.fail = !1, o.zoom = b.zoom, o.center = b.center, o.active_alternative = b.active_alternative, o.name = b.destination_name;
                    for (var p = 0; p < m.length; p++) o.positions.push(new L.LatLng(0, 0));
                    for (var p = 0; p < m.length; p++) {
                        var q = m[p];
                        if (q.match(/^\s*[-+]?[0-9]*\.?[0-9]+\s*[,;]\s*[-+]?[0-9]*\.?[0-9]+\s*$/)) {
                            var r = q.split(/[,;]/);
                            OSRM.Geocoder._showInitResults([{
                                lat: r[0],
                                lon: r[1]
                            }], {
                                id: p,
                                callback: n
                            })
                        } else {
                            OSRM.GUI.exclusiveNotify(OSRM.loc("NOTIFICATION_GEOCODERWAIT_HEADER"), OSRM.loc("NOTIFICATION_GEOCODERWAIT_BODY"), !1);
                            var s = OSRM.DEFAULTS.HOST_GEOCODER_URL + "?format=json&json_callback=%jsonp" + OSRM.DEFAULTS.GEOCODER_BOUNDS + "&accept-language=" + OSRM.Localization.current_language + "&limit=1&q=" + q;
                            OSRM.JSONP.call(s, OSRM.Geocoder._showInitResults, OSRM.Geocoder._showInitResults, OSRM.DEFAULTS.JSONP_TIMEOUT, "init_geocoder_" + p, {
                                id: p,
                                callback: n
                            })
                        }
                    }
                } else;
        }
    }, OSRM.Browser.onLoadHandler(OSRM.init), OSRM.Control.Attribution = L.Control.extend({
        options: {
            position: "bottomright",
            prefix: 'Powered by <a href="http://leaflet.cloudmade.com">Leaflet</a>',
            postfix: ""
        },
        initialize: function(a) {
            L.Util.setOptions(this, a), this._attributions = {}
        },
        onAdd: function(a) {
            return this._container = L.DomUtil.create("div", "leaflet-control-attribution"), L.DomEvent.disableClickPropagation(this._container), a.on("layeradd", this._onLayerAdd, this).on("layerremove", this._onLayerRemove, this), this._update(), this._container
        },
        onRemove: function(a) {
            a.off("layeradd", this._onLayerAdd).off("layerremove", this._onLayerRemove)
        },
        setPrefix: function(a) {
            return this.options.prefix = a, this._update(), this
        },
        setPostfix: function(a) {
            return this.options.postfix = a, this._update(), this
        },
        addAttribution: function(a) {
            return a ? (this._attributions[a] || (this._attributions[a] = 0), this._attributions[a] ++, this._update(), this) : void 0
        },
        removeAttribution: function(a) {
            return a ? (this._attributions[a] --, this._update(), this) : void 0
        },
        _update: function() {
            if (this._map) {
                var a = [];
                for (var b in this._attributions) this._attributions.hasOwnProperty(b) && this._attributions[b] && a.push(b);
                var c = [];
                this.options.prefix && c.push(this.options.prefix), a.length && c.push(a.join(", ")), this.options.postfix && c.push(this.options.postfix), this._container.innerHTML = c.join(" &#8212; ")
            }
        },
        _onLayerAdd: function(a) {
            a.layer.getAttribution && this.addAttribution(a.layer.getAttribution())
        },
        _onLayerRemove: function(a) {
            a.layer.getAttribution && this.removeAttribution(a.layer.getAttribution())
        }
    }), OSRM.Control.Layers = L.Control.Layers.extend({
        getActiveLayerName: function() {
            var a, b, c, d = this._form.getElementsByTagName("input"),
                e = d.length;
            for (a = 0; e > a; a++)
                if (b = d[a], c = this._layers[b.layerId], b.checked && !c.overlay) return c.name
        },
        getActiveLayer: function() {
            var a, b, c, d = this._form.getElementsByTagName("input"),
                e = d.length;
            for (a = 0; e > a; a++)
                if (b = d[a], c = this._layers[b.layerId], b.checked && !c.overlay) return c.layer
        },
        setActiveLayerByName: function(a) {
            var b, c, d, e = this._form.getElementsByTagName("input"),
                f = e.length;
            for (b = 0; f > b; b++) c = e[b], d = this._layers[c.layerId], c.checked = d.name != a ? !1 : !0;
            this._onInputClick()
        },
        setLayerLabels: function() {
            var a, b, c = this._form.getElementsByTagName("input"),
                d = c.length;
            for (tileServers = OSRM.DEFAULTS.TILE_SERVERS.length, a = 0; d > a; a++)
                if (b = c[a], tileServers > a) b.parentNode.lastChild.textContent = OSRM.loc("TILE_SERVER_" + a) == "TILE_SERVER_" + a ? " " + OSRM.DEFAULTS.TILE_SERVERS[a].display_name : " " + OSRM.loc("TILE_SERVER_" + a);
                else {
                    var e = a - tileServers;
                    b.parentNode.lastChild.textContent = OSRM.loc("OVERLAY_SERVER_" + e) == "OVERLAY_SERVER_" + e ? " " + OSRM.DEFAULTS.OVERLAY_SERVERS[e].display_name : " " + OSRM.loc("OVERLAY_SERVER_" + e)
                }
        },
        _initLayout: function() {
            L.Control.Layers.prototype._initLayout.apply(this), this._container.className = "box-wrapper gui-control-wrapper", this._layersLink.className = "box-content gui-control gui-layers", this._form.className = "box-content gui-control gui-layers-list medium-font", this._baseLayersList.className = "gui-layers-base", this._separator.className = "gui-layers-separator", this._overlaysList.className = "gui-layers-overlays"
        },
        _expand: function() {
            L.DomUtil.addClass(this._container, "gui-layers-expanded")
        },
        _collapse: function() {
            this._container.className = this._container.className.replace(" gui-layers-expanded", "")
        }
    }), OSRM.Control.Locations = L.Control.extend({
        options: {
            position: "topright"
        },
        onAdd: function(a) {
            var b = L.DomUtil.create("div", "box-wrapper gui-control-wrapper");
            return L.DomEvent.disableClickPropagation(b), this._userButton = this._createButton("gui-locations-user", b, OSRM.GUI.zoomOnUser, a, !!navigator.geolocation), this._routeButton = this._createButton("gui-locations-route", b, OSRM.GUI.zoomOnRoute, a, !1), this._container = b, b
        },
        _createButton: function(a, b, c, d, e) {
            var f = 0 == e ? "-inactive" : "",
                g = "box-content gui-control" + f + " " + a + f,
                h = L.DomUtil.create("a", g, b);
            return h.title = a, L.DomEvent.on(h, "click", L.DomEvent.stopPropagation).on(h, "click", L.DomEvent.preventDefault).on(h, "click", c, d).on(h, "dblclick", L.DomEvent.stopPropagation), h
        },
        activateRoute: function() {
            this._routeButton.className = "box-content gui-control gui-locations-route"
        },
        deactivateRoute: function() {
            this._routeButton.className = "box-content gui-control-inactive gui-locations-route-inactive"
        },
        setTooltips: function(a, b) {
            this._userButton.title = a, this._routeButton.title = b
        }
    }), OSRM.Control.Zoom = L.Control.extend({
        options: {
            position: "topleft"
        },
        onAdd: function(a) {
            var b = "gui-zoom",
                c = L.DomUtil.create("div", "box-wrapper gui-control-wrapper");
            return L.DomEvent.disableClickPropagation(c), this._map = a, this._zoomInButton = this._createButton("", "Zoom in", b + "-in", c, this._zoomIn, this), this._zoomOutButton = this._createButton("", "Zoom out", b + "-out", c, this._zoomOut, this), this._container = c, a.on("zoomend zoomlevelschange", this._updateDisabled, this), c
        },
        onRemove: function(a) {
            a.off("zoomend zoomlevelschange", this._updateDisabled, this)
        },
        _zoomIn: function(a) {
            this._map.zoomIn(a.shiftKey ? 3 : 1)
        },
        _zoomOut: function(a) {
            this._map.zoomOut(a.shiftKey ? 3 : 1)
        },
        _createButton: function(a, b, c, d, e, f) {
            var g = "box-content gui-control " + c,
                h = L.DomUtil.create("a", g, d);
            h.innerHTML = a, h.title = b;
            var i = L.DomEvent.stopPropagation;
            return L.DomEvent.on(h, "click", i).on(h, "mousedown", i).on(h, "dblclick", i).on(h, "click", L.DomEvent.preventDefault).on(h, "click", e, f), h
        },
        _updateDisabled: function() {
            var a = this._map,
                b = "gui-control",
                c = "gui-zoom";
            a._zoom === a.getMinZoom() ? (L.DomUtil.changeClass(this._zoomOutButton, b, b + "-inactive"), L.DomUtil.changeClass(this._zoomOutButton, c + "-out", c + "-out-inactive")) : (L.DomUtil.changeClass(this._zoomOutButton, b + "-inactive", b), L.DomUtil.changeClass(this._zoomOutButton, c + "-out-inactive", c + "-out")), a._zoom === a.getMaxZoom() ? (L.DomUtil.changeClass(this._zoomInButton, b, b + "-inactive"), L.DomUtil.changeClass(this._zoomInButton, c + "-in", c + "-in-inactive")) : (L.DomUtil.changeClass(this._zoomInButton, b + "-inactive", b), L.DomUtil.changeClass(this._zoomInButton, c + "-in-inactive", c + "-in"))
        },
        hide: function() {
            this._container && (this._container.style.visibility = "hidden")
        },
        hide: function() {
            this._container && (this._container.style.visibility = "hidden")
        },
        show: function() {
            this._container && (this._container.style.top = "5px", this._container.style.left = (1 == OSRM.G.main_handle.boxVisible() ? OSRM.G.main_handle.boxWidth() + 10 : "30") + "px", this._container.style.visibility = "visible")
        },
        setTooltips: function(a, b) {
            this._zoomInButton.title = a, this._zoomOutButton.title = b
        }
    }), OSRM.Control.Map = L.Map.extend({
        _boundsInsideView: function(a) {
            var b = this.getBounds(),
                c = this.project(b.getSouthWest()),
                d = this.project(b.getNorthEast()),
                e = this.project(a.getSouthWest()),
                f = this.project(a.getNorthEast());
            return d.y > f.y ? !1 : d.x < f.x ? !1 : c.y < e.y ? !1 : c.x > e.x ? !1 : !0
        },
        setViewBounds: function(a) {
            var b = this.getBoundsZoom(a);
            this._zoom > b ? this.setView(a.getCenter(), b) : this._boundsInsideView(a) || this.setView(a.getCenter(), this._zoom)
        },
        setViewUI: function(a, b, c) {
            if (OSRM.G.main_handle.boxVisible()) {
                var d = this.project(a, b);
                d.x -= OSRM.G.main_handle.boxWidth() / 2, a = this.unproject(d, b)
            }
            this.setView(a, b, c)
        },
        setViewBoundsUI: function(a) {
            var b = a.getSouthWest(),
                c = a.getNorthEast(),
                d = this.getBoundsZoom(a),
                e = this.project(b, d);
            e.x -= OSRM.G.main_handle.boxVisible() ? OSRM.G.main_handle.boxWidth() + 20 : 20, e.y += 20;
            var f = this.project(c, d);
            f.y -= 20, f.x += 20, a.extend(this.unproject(e, d)), a.extend(this.unproject(f, d)), this.setViewBounds(a)
        },
        fitBoundsUI: function(a) {
            var b = a.getSouthWest(),
                c = a.getNorthEast(),
                d = this.getBoundsZoom(a),
                e = this.project(b, d);
            e.x -= OSRM.G.main_handle.boxVisible() ? OSRM.G.main_handle.boxWidth() + 20 : 20, e.y += 20;
            var f = this.project(c, d);
            f.y -= 20, f.x += 20, a.extend(this.unproject(e, d)), a.extend(this.unproject(f, d)), this.fitBounds(a)
        },
        getBoundsUI: function() {
            var a = this.getPixelBounds();
            OSRM.G.main_handle.boxVisible() && (a.min.x += OSRM.G.main_handle.boxWidth());
            var b = this.unproject(new L.Point(a.min.x, a.max.y), this._zoom, !0),
                c = this.unproject(new L.Point(a.max.x, a.min.y), this._zoom, !0);
            return new L.LatLngBounds(b, c)
        },
        getCenterUI: function(a) {
            var b = this.getSize();
            OSRM.G.main_handle.boxVisible() && (b.x += OSRM.G.main_handle.boxWidth());
            var c = this._getTopLeftPoint().add(b.divideBy(2));
            return this.unproject(c, this._zoom, a)
        },
        getActiveLayerId: function() {
            for (var a = 0, b = OSRM.DEFAULTS.TILE_SERVERS, c = this.layerControl.getActiveLayerName(), d = 0, e = b.length; e > d; d++)
                if (b[d].display_name == c) {
                    a = d;
                    break
                }
            return a
        }
    }), OSRM.Marker = function(a, b, c) {
        this.label = a ? a : "marker", this.position = c ? c : new L.LatLng(0, 0), this.description = null, this.marker = new L.LabelMarker(this.position, b), this.marker.parent = this, this.shown = !1, this.hint = null
    }, OSRM.extend(OSRM.Marker, {
        show: function() {
            OSRM.G.map.addLayer(this.marker), this.shown = !0
        },
        hide: function() {
            if (OSRM.G.map.removeLayer(this.marker), this.shown = !1, "highlight" == this.label && this.description) {
                var a = document.getElementById("description-" + this.description);
                a && (a.className = "description-body-item"), this.description = null
            }
        },
        setPosition: function(a) {
            this.position = a, this.marker.setLatLng(a), this.hint = null
        },
        getPosition: function() {
            return this.position
        },
        getLat: function() {
            return this.position.lat
        },
        getLng: function() {
            var a = this.position.lng - 360 * Math.floor(this.position.lng / 360 + .5);
            return a
        },
        isShown: function() {
            return this.shown
        },
        centerView: function(a) {
            void 0 == a && (a = OSRM.DEFAULTS.ZOOM_LEVEL), OSRM.G.map.setViewUI(this.position, a)
        },
        toString: function() {
            return 'OSRM.Marker: "' + this.label + '", ' + this.position + ")"
        }
    }), OSRM.RouteMarker = function(a, b) {
        b.baseicon = b.icon, OSRM.RouteMarker.prototype.base.constructor.apply(this, arguments), this.label = a ? a : "route_marker", this.marker.on("click", this.onClick), this.marker.on("drag", this.onDrag), this.marker.on("dragstart", this.onDragStart), this.marker.on("dragend", this.onDragEnd)
    }, OSRM.inheritFrom(OSRM.RouteMarker, OSRM.Marker), OSRM.extend(OSRM.RouteMarker, {
        onClick: function(a) {
            if (1 != a.originalEvent.shiftKey && 1 != a.originalEvent.metaKey && 1 != a.originalEvent.altKey) {
                for (var b = 0; b < OSRM.G.markers.route.length; b++)
                    if (OSRM.G.markers.route[b].marker === this) {
                        OSRM.G.markers.removeMarker(b);
                        break
                    }
                OSRM.Routing.getRoute(), OSRM.G.markers.highlight.hide(), OSRM.G.markers.dragger.hide()
            }
        },
        onDrag: function(a) {
            this.parent.setPosition(a.target.getLatLng()), OSRM.G.markers.route.length > 1 && OSRM.Routing.getRoute_Dragging(), OSRM.Geocoder.updateLocation(this.parent.label)
        },
        onDragStart: function() {
            OSRM.GUI.deactivateTooltip("DRAGGING"), OSRM.G.dragging = !0, this.changeIcon(this.options.dragicon), this.parent.description = null;
            for (var a = 0; a < OSRM.G.markers.route.length; a++)
                if (OSRM.G.markers.route[a].marker === this) {
                    OSRM.G.dragid = a;
                    break
                }
            this.parent != OSRM.G.markers.highlight && OSRM.G.markers.highlight.hide(), this.parent != OSRM.G.markers.dragger && OSRM.G.markers.dragger.hide(), OSRM.G.route.isShown() && OSRM.G.route.showOldRoute()
        },
        onDragEnd: function(a) {
            OSRM.G.dragging = !1, this.changeIcon(this.options.baseicon), this.parent.setPosition(a.target.getLatLng()), OSRM.G.route.isShown() ? (OSRM.Routing.getRoute(), OSRM.G.route.hideOldRoute(), OSRM.G.route.hideUnnamedRoute()) : (OSRM.Geocoder.updateAddress(this.parent.label), OSRM.GUI.clearResults())
        },
        toString: function() {
            return 'OSRM.RouteMarker: "' + this.label + '", ' + this.position + ")"
        }
    }), OSRM.DragMarker = function(a) {
        OSRM.DragMarker.prototype.base.constructor.apply(this, arguments), this.label = a ? a : "drag_marker"
    }, OSRM.inheritFrom(OSRM.DragMarker, OSRM.RouteMarker), OSRM.extend(OSRM.DragMarker, {
        onClick: function(a) {
            if (this.parent != OSRM.G.markers.dragger) this.parent.hide();
            else {
                var b = OSRM.Via.findViaIndex(a.target.getLatLng());
                OSRM.G.markers.route.splice(b + 1, 0, this.parent), OSRM.RouteMarker.prototype.onDragStart.call(this, a), OSRM.G.markers.route[OSRM.G.dragid] = new OSRM.RouteMarker(OSRM.C.VIA_LABEL, {
                    draggable: !0,
                    icon: OSRM.G.icons["marker-via"],
                    dragicon: OSRM.G.icons["marker-via-drag"]
                }, a.target.getLatLng()), OSRM.G.markers.route[OSRM.G.dragid].show(), OSRM.RouteMarker.prototype.onDragEnd.call(this, a), this.parent.hide()
            }
        },
        onDragStart: function(a) {
            var b = OSRM.Via.findViaIndex(a.target.getLatLng());
            OSRM.G.markers.route.splice(b + 1, 0, this.parent), OSRM.RouteMarker.prototype.onDragStart.call(this, a)
        },
        onDragEnd: function(a) {
            OSRM.G.markers.route[OSRM.G.dragid] = new OSRM.RouteMarker(OSRM.C.VIA_LABEL, {
                draggable: !0,
                icon: OSRM.G.icons["marker-via"],
                dragicon: OSRM.G.icons["marker-via-drag"]
            }, a.target.getLatLng()), OSRM.G.markers.route[OSRM.G.dragid].show(), OSRM.RouteMarker.prototype.onDragEnd.call(this, a), this.parent.hide()
        },
        toString: function() {
            return 'OSRM.DragMarker: "' + this.label + '", ' + this.position + ")"
        }
    }), OSRM.SimpleRoute = function(a, b) {
        this.label = a ? a : "route", this.route = new L.Polyline([], b), this.shown = !1, this.route.on("click", this.onClick)
    }, OSRM.extend(OSRM.SimpleRoute, {
        show: function() {
            OSRM.G.map.addLayer(this.route), this.shown = !0
        },
        hide: function() {
            OSRM.G.map.removeLayer(this.route), this.shown = !1
        },
        isShown: function() {
            return this.shown
        },
        getPoints: function() {
            return this.route._originalPoints
        },
        getPositions: function() {
            return this.route.getLatLngs()
        },
        setPositions: function(a) {
            this.route.setLatLngs(a)
        },
        setStyle: function(a) {
            this.route.setStyle(a)
        },
        centerView: function() {
            var a = new L.LatLngBounds(this.getPositions());
            OSRM.g.map.fitBoundsUI(a)
        },
        onClick: function(a) {
            if (1 != a.originalEvent.shiftKey && 1 != a.originalEvent.metaKey && 1 != a.originalEvent.altKey) {
                var b = Math.max(0, OSRM.Via.findViaIndex(a.latlng)),
                    c = OSRM.G.markers.setVia(b, a.latlng);
                OSRM.G.markers.route[c].show(), OSRM.Routing.getRoute()
            }
        },
        toString: function() {
            return "OSRM.Route(" + this.label + ", " + this.route.getLatLngs().length + " points)"
        }
    }), OSRM.MultiRoute = function(a) {
        this.label = a ? a : "multiroute", this.route = new L.LayerGroup, this.shown = !1
    }, OSRM.extend(OSRM.MultiRoute, {
        show: function() {
            OSRM.G.map.addLayer(this.route), this.shown = !0
        },
        hide: function() {
            OSRM.G.map.removeLayer(this.route), this.shown = !1
        },
        isShown: function() {
            return this.shown
        },
        addRoute: function(a) {
            var b = new L.Polyline(a);
            b.on("click", function(a) {
                OSRM.G.route.fire("click", a)
            }), this.route.addLayer(b)
        },
        clearRoutes: function() {
            this.route.clearLayers()
        },
        setStyle: function(a) {
            this.route.invoke("setStyle", a)
        },
        toString: function() {
            return "OSRM.MultiRoute(" + this.label + ")"
        }
    }), OSRM.GLOBALS.map = null, OSRM.GLOBALS.localizable_maps = [], OSRM.Map = {
        init: function() {
            null == OSRM.G.main_handle && OSRM.GUI.init();
            for (var a = OSRM.DEFAULTS.TILE_SERVERS, b = {}, c = 0, d = a.length; d > c; c++) a[c].options.attribution = a[c].attribution, b[a[c].display_name] = new L.TileLayer(a[c].url, a[c].options), L.Util.stamp(b[a[c].display_name]);
            for (var e = OSRM.DEFAULTS.OVERLAY_SERVERS, f = {}, c = 0, d = e.length; d > c; c++) e[c].options.attribution = e[c].attribution, f[e[c].display_name] = new L.TileLayer(e[c].url, e[c].options), L.Util.stamp(f[e[c].display_name]);
            OSRM.G.map = new OSRM.Control.Map("map", {
                center: new L.LatLng(OSRM.DEFAULTS.ONLOAD_LATITUDE, OSRM.DEFAULTS.ONLOAD_LONGITUDE),
                zoom: OSRM.DEFAULTS.ONLOAD_ZOOM_LEVEL,
                layers: [],
                zoomAnimation: !1,
                fadeAnimation: !1,
                zoomControl: !1,
                attributionControl: !1,
                worldCopyJump: !0
            }), OSRM.G.map.attributionControl = new OSRM.Control.Attribution, OSRM.G.map.attributionControl.addTo(OSRM.G.map), OSRM.G.map.locationsControl = new OSRM.Control.Locations, OSRM.G.map.locationsControl.addTo(OSRM.G.map), OSRM.G.map.addLayer(b[a[0].display_name]), OSRM.G.map.layerControl = new OSRM.Control.Layers(b, f), OSRM.G.map.layerControl.addTo(OSRM.G.map), OSRM.G.map.zoomControl = new OSRM.Control.Zoom, OSRM.G.map.zoomControl.addTo(OSRM.G.map), OSRM.G.map.zoomControl.show(), OSRM.G.map.scaleControl = new L.Control.Scale, OSRM.G.map.scaleControl.options.metric = 1 != OSRM.G.DISTANCE_FORMAT, OSRM.G.map.scaleControl.options.imperial = 1 == OSRM.G.DISTANCE_FORMAT, OSRM.G.map.scaleControl.addTo(OSRM.G.map), OSRM.G.map.on("zoomend", OSRM.Map.zoomed), OSRM.G.map.on("click", OSRM.Map.click), OSRM.G.map.on("contextmenu", OSRM.Map.contextmenu), OSRM.G.map.on("mousemove", OSRM.Map.mousemove)
        },
        initFinally: function() {
            L.Util.setOptions(OSRM.G.map, {
                zoomAnimation: !0,
                fadeAnimation: !0
            })
        },
        initPosition: function() {
            var a = new L.LatLng(OSRM.DEFAULTS.ONLOAD_LATITUDE, OSRM.DEFAULTS.ONLOAD_LONGITUDE);
            OSRM.G.map.setViewUI(a, OSRM.DEFAULTS.ONLOAD_ZOOM_LEVEL, !0), navigator.geolocation && -1 == document.URL.indexOf("file://") && navigator.geolocation.getCurrentPosition(OSRM.Map.geolocationResponse)
        },
        zoomed: function() {
            var a = OSRM.G.route.getZoomLevel() - OSRM.G.map.getZoom();
            a >= 0 && 3 >= a || (OSRM.G.dragging ? OSRM.Routing.getRoute_Dragging() : OSRM.Routing.getRoute_Redraw({
                keepAlternative: !0
            }))
        },
        contextmenu: function() {},
        mousemove: function(a) {
            //OSRM.Via.drawDragMarker(a)
        },
        click: function(e) {
			OSRM.GUI.deactivateTooltip( "CLICKING" );
			if( e.originalEvent.shiftKey==true || e.originalEvent.metaKey==true || e.originalEvent.altKey==true )	// only create/remove markers on simple clicks
				return;
			if( !OSRM.G.markers.hasSource() ) {
				var index = OSRM.G.markers.setSource( e.latlng );
				OSRM.Geocoder.updateAddress( OSRM.C.SOURCE_LABEL, OSRM.C.DO_FALLBACK_TO_LAT_LNG );
				OSRM.G.markers.route[index].show();		
				OSRM.Routing.getRoute( {recenter:OSRM.G.markers.route.length == 2} );	// allow recentering when the route is first shown 
			} else if( !OSRM.G.markers.hasTarget() ) {
				var index = OSRM.G.markers.setTarget( e.latlng );
				OSRM.Geocoder.updateAddress( OSRM.C.TARGET_LABEL, OSRM.C.DO_FALLBACK_TO_LAT_LNG );
				OSRM.G.markers.route[index].show();
				OSRM.Routing.getRoute( {recenter:OSRM.G.markers.route.length == 2} );	// allow recentering when the route is first shown
			} else {
				var index = OSRM.G.markers.setVia(OSRM.G.markers.route.length - 2, e.latlng );
				OSRM.G.markers.route[index].show();
				OSRM.Routing.getRoute( {recenter:false} );
			}
        },
        geolocationResponse: function(a) {
            var b = new L.LatLng(a.coords.latitude, a.coords.longitude);
            OSRM.G.map.setViewUI(b, OSRM.DEFAULTS.ZOOM_LEVEL)
        }
    }, OSRM.Markers = function() {
        this.route = new Array, this.highlight = new OSRM.DragMarker("highlight", {
            zIndexOffset: -1,
            draggable: !0,
            icon: OSRM.G.icons["marker-highlight"],
            dragicon: OSRM.G.icons["marker-highlight-drag"]
        }), this.hover = new OSRM.Marker("hover", {
            zIndexOffset: -1,
            draggable: !1,
            icon: OSRM.G.icons["marker-highlight"]
        }), this.dragger = new OSRM.DragMarker("drag", {
            draggable: !0,
            icon: OSRM.G.icons["marker-drag"],
            dragicon: OSRM.G.icons["marker-drag"]
        }), this.initial_vias = new Array
    }, OSRM.extend(OSRM.Markers, {
        reset: function() {
            for (var a = 0; a < this.route.length; a++) this.route[a].hide();
            this.route.splice(0, this.route.length), document.getElementById("gui-delete-source").style.visibility = "hidden", document.getElementById("gui-delete-target").style.visibility = "hidden", this.highlight.hide(), this.dragger.hide(), this.initial_vias.length = 0
        },
        removeVias: function() {
            for (var a = 1; a < this.route.length - 1; a++) this.route[a].hide();
            this.route.splice(1, this.route.length - 2)
        },
        setSource: function(a) {
            if (this.route[0] && this.route[0].label == OSRM.C.SOURCE_LABEL ? this.route[0].setPosition(a) : this.route.splice(0, 0, new OSRM.RouteMarker(OSRM.C.SOURCE_LABEL, {
                    draggable: !0,
                    icon: OSRM.G.icons["marker-source"],
                    dragicon: OSRM.G.icons["marker-source-drag"]
                }, a)), document.getElementById("gui-delete-source").style.visibility = "visible", this.initial_vias.length > 0) {
                for (var b = 0; b < this.initial_vias.length; b++) OSRM.G.markers.setVia(b, this.initial_vias[b]);
                for (var b = 1; b < OSRM.G.markers.route.length - 1; b++) OSRM.G.markers.route[b].show()
            }
            return 0
        },
        setTarget: function(a) {
            if (this.route[this.route.length - 1] && this.route[this.route.length - 1].label == OSRM.C.TARGET_LABEL ? this.route[this.route.length - 1].setPosition(a) : this.route.splice(this.route.length, 0, new OSRM.RouteMarker(OSRM.C.TARGET_LABEL, {
                    draggable: !0,
                    icon: OSRM.G.icons["marker-target"],
                    dragicon: OSRM.G.icons["marker-target-drag"]
                }, a)), document.getElementById("gui-delete-target").style.visibility = "visible", this.initial_vias.length > 0) {
                for (var b = 0; b < this.initial_vias.length; b++) OSRM.G.markers.setVia(b, this.initial_vias[b]);
                for (var b = 1; b < OSRM.G.markers.route.length - 1; b++) OSRM.G.markers.route[b].show()
            }
            return this.route.length - 1
        },
        setVia: function(a, b) {
            return this.route.length < 2 || a > this.route.length - 2 ? -1 : (this.route.splice(a + 1, 0, new OSRM.RouteMarker(OSRM.C.VIA_LABEL, {
                draggable: !0,
                icon: OSRM.G.icons["marker-via"],
                dragicon: OSRM.G.icons["marker-via-drag"]
            }, b)), a + 1)
        },
        removeMarker: function(a) {
            a >= this.route.length || (0 == a && this.route[0].label == OSRM.C.SOURCE_LABEL ? (this.removeVias(), document.getElementById("gui-input-source").value = "", document.getElementById("gui-delete-source").style.visibility = "hidden", OSRM.GUI.clearResults()) : a == this.route.length - 1 && this.route[this.route.length - 1].label == OSRM.C.TARGET_LABEL && (this.removeVias(), a = this.route.length - 1, document.getElementById("gui-input-target").value = "", document.getElementById("gui-delete-target").style.visibility = "hidden", OSRM.GUI.clearResults()), this.route[a].hide(), this.route.splice(a, 1), this.initial_vias.length = 0)
        },
        reverseDescriptions: function() {
            for (var a = this.route.length - 1, b = this.route.length / 2, c = 0; b > c; ++c) {
                var d = this.route[c].description;
                this.route[c].description = this.route[a - c].description, this.route[a - c].description = d
            }
        },
        reverseMarkers: function() {
            var a = this.route.length,
                b = this.route[0].getPosition();
            this.route[0].setPosition(this.route[a - 1].getPosition()), this.route[a - 1].setPosition(b);
            var c = this.route[0];
            this.route[0] = this.route[a - 1], this.route[a - 1] = c, this.route.reverse(), OSRM.GUI.clearResults(), this.reverseInitialVias()
        },
        hasSource: function() {
            return this.route[0] && this.route[0].label == OSRM.C.SOURCE_LABEL ? !0 : !1
        },
        hasTarget: function() {
            return this.route[this.route.length - 1] && this.route[this.route.length - 1].label == OSRM.C.TARGET_LABEL ? !0 : !1
        },
        addInitialVia: function(a) {
            this.initial_vias.push(a)
        },
        reverseInitialVias: function() {
            this.initial_vias.reverse()
        },
        relabelViaMarkers: function() {
            for (var a = 0, b = this.route.length; b > a; a++) this.route[a].marker.setLabel(a)
        }
    }), OSRM.Route = function() {
        this._current_route = new OSRM.SimpleRoute("current", {
            dashArray: ""
        }), this._alternative_route = new OSRM.SimpleRoute("alternative", {
            dashArray: ""
        }), this._old_route = new OSRM.SimpleRoute("old", {
            color: "#123",
            dashArray: ""
        }), this._unnamed_route = new OSRM.MultiRoute("unnamed"), this._current_route_style = {
            color: "#0033FF",
            weight: 5,
            dashArray: ""
        }, this._current_noroute_style = {
            color: "#222222",
            weight: 2,
            dashArray: "8,6"
        }, this._old_route_style = {
            color: "#112233",
            weight: 5,
            dashArray: ""
        }, this._old_noroute_style = {
            color: "#000000",
            weight: 2,
            dashArray: "8,6"
        }, this._unnamed_route_style = {
            color: "#FF00FF",
            weight: 10,
            dashArray: ""
        }, this._old_unnamed_route_style = {
            color: "#990099",
            weight: 10,
            dashArray: ""
        }, this._alternative_route_style = {
            color: "#770033",
            weight: 5,
            opacity: .6,
            dashArray: ""
        }, this._noroute = OSRM.Route.ROUTE, this._history = new OSRM.HistoryRoute, this._zoomlevel = 0
    }, OSRM.Route.NOROUTE = !0, OSRM.Route.ROUTE = !1, OSRM.extend(OSRM.Route, {
        showRoute: function(a, b) {
            this._noroute = b, this._current_route.setPositions(a), this._current_route.setStyle(this._noroute == OSRM.Route.NOROUTE ? this._current_noroute_style : this._current_route_style), this._current_route.show(), this._history.fetchHistoryRoute(), this._history.showHistoryRoutes(), this._history.storeHistoryRoute(), this._zoomlevel = OSRM.G.map.getZoom()
        },
        hideRoute: function() {
            this._current_route.hide(), this._unnamed_route.hide(), this._history.fetchHistoryRoute(), this._history.showHistoryRoutes(), this._zoomlevel = 0, OSRM.GUI.deactivateRouteFeatures()
        },
        showUnnamedRoute: function(a) {
            this._unnamed_route.clearRoutes();
            for (var b = 0; b < a.length; b++) this._unnamed_route.addRoute(a[b]);
            this._unnamed_route.setStyle(this._unnamed_route_style), this._unnamed_route.show()
        },
        hideUnnamedRoute: function() {
            this._unnamed_route.hide()
        },
        _raiseUnnamedRoute: function() {
            this._unnamed_route.isShown() && (this._unnamed_route.hide(), this._unnamed_route.show())
        },
        showOldRoute: function() {
            this._old_route.setPositions(this._current_route.getPositions()), this._old_route.setStyle(this._noroute == OSRM.Route.NOROUTE ? this._old_noroute_style : this._old_route_style), this._old_route.show(), this._raiseUnnamedRoute(), this._unnamed_route.setStyle(this._old_unnamed_route_style)
        },
        hideOldRoute: function() {
            this._old_route.hide()
        },
        showAlternativeRoute: function(a) {
            this._alternative_route.setPositions(a), this._alternative_route.setStyle(this._alternative_route_style), this._alternative_route.show()
        },
        hideAlternativeRoute: function() {
            this._alternative_route.hide()
        },
        isShown: function() {
            return this._current_route.isShown()
        },
        isRoute: function() {
            return !this._noroute
        },
        getPositions: function() {
            return this._current_route.getPositions()
        },
        getPoints: function() {
            return this._current_route.getPoints()
        },
        getZoomLevel: function() {
            return this._zoomlevel
        },
        reset: function() {
            this.hideRoute(), this._old_route.hide(), this._noroute = OSRM.Route.ROUTE, this._history.clearHistoryRoutes()
        },
        fire: function(a, b) {
            this._current_route.route.fire(a, b)
        },
        centerView: function() {
            this._current_route.centerView()
        },
        activateHistoryRoutes: function() {
            this._history.activate()
        },
        deactivateHistoryRoutes: function() {
            this._history.deactivate()
        }
    }), OSRM.HistoryRoute = function() {
        this._history_styles = [{
            color: "#FFFFFF",
            opacity: .5,
            weight: 5,
            dashArray: ""
        }, {
            color: "#0000DD",
            opacity: .45,
            weight: 5,
            dashArray: ""
        }, {
            color: "#0000BB",
            opacity: .4,
            weight: 5,
            dashArray: ""
        }, {
            color: "#000099",
            opacity: .35,
            weight: 5,
            dashArray: ""
        }, {
            color: "#000077",
            opacity: .3,
            weight: 5,
            dashArray: ""
        }, {
            color: "#000055",
            opacity: .25,
            weight: 5,
            dashArray: ""
        }, {
            color: "#000033",
            opacity: .2,
            weight: 5,
            dashArray: ""
        }, {
            color: "#000011",
            opacity: .15,
            weight: 5,
            dashArray: ""
        }, {
            color: "#000000",
            opacity: .1,
            weight: 5,
            dashArray: ""
        }], this._history_length = this._history_styles.length, this._history = [];
        for (var a = 0, b = this._history_length; b > a; a++) {
            var c = {};
            c.route = new OSRM.SimpleRoute("current", {
                dashArray: ""
            }), c.markers = [], c.checksum = null, this._history.push(c)
        }
        this._initiate_redrawHistory = OSRM.bind(this, this._getRoute_RedrawHistory), this._callback_redrawHistory = OSRM.bind(this, this._showRoute_RedrawHistory)
    }, OSRM.extend(OSRM.HistoryRoute, {
        activate: function() {
            this.storeHistoryRoute = this._storeHistoryRoute, this.fetchHistoryRoute = this._fetchHistoryRoute, this.showHistoryRoutes = this._showHistoryRoutes, this.clearHistoryRoutes = this._clearHistoryRoutes, OSRM.G.map.on("zoomend", this._initiate_redrawHistory), this.storeHistoryRoute()
        },
        deactivate: function() {
            this.clearHistoryRoutes(), this.storeHistoryRoute = this.empty, this.fetchHistoryRoute = this.empty, this.showHistoryRoutes = this.empty, this.clearHistoryRoutes = this.empty, OSRM.G.map.off("zoomend", this._initiate_redrawHistory)
        },
        empty: function() {},
        storeHistoryRoute: function() {},
        fetchHistoryRoute: function() {},
        showHistoryRoutes: function() {},
        clearHistoryRoutes: function() {},
        _storeHistoryRoute: function() {
            var a = OSRM.G.route;
            if (a.isShown() && a.isRoute()) {
                var b = OSRM.G.response.hint_data;
                this._history[0].route.setPositions(a.getPositions()), this._history[0].checksum = b.checksum, this._history[0].markers = [];
                for (var c = this._getCurrentMarkers(), d = 0, e = c.length; e > d; d++) {
                    var f = {
                        lat: c[d].lat,
                        lng: c[d].lng,
                        hint: b.locations[d]
                    };
                    this._history[0].markers.push(f)
                }
            }
        },
        _fetchHistoryRoute: function() {
            if (0 != this._history[0].markers.length && !(OSRM.G.route.isShown() && this._equalMarkers(this._history[0].markers, this._getCurrentMarkers()) || this._equalMarkers(this._history[0].markers, this._history[1].markers))) {
                for (var a = this._history_length - 1; a > 0; a--) this._history[a].route.setPositions(this._history[a - 1].route.getPositions()), this._history[a].markers = this._history[a - 1].markers, this._history[a].checksum = this._history[a - 1].checksum;
                this._history[0].route.setPositions([]), this._history[0].markers = [], this._history[0].checksum = null
            }
        },
        _showHistoryRoutes: function() {
            for (var a = 1, b = this._history_length; b > a; a++) this._history[a].route.setStyle(this._history_styles[a]), this._history[a].route.show(), OSRM.G.route.hideOldRoute()
        },
        _clearHistoryRoutes: function() {
            for (var a = 0, b = this._history_length; b > a; a++) this._history[a].route.hide(), this._history[a].route.setPositions([]), this._history[a].markers = [], this._history[a].checksum = null
        },
        _getCurrentMarkers: function() {
            var a = [],
                b = OSRM.G.route.getPositions();
            if (0 == b.length) return a;
            for (var c = 0; c < OSRM.G.response.via_points.length; c++) a.push({
                lat: OSRM.G.response.via_points[c][0],
                lng: OSRM.G.response.via_points[c][1]
            });
            return a
        },
        _equalMarkers: function(a, b) {
            if (a.length != b.length) return !1;
            for (var c = OSRM.C.PRECISION, d = 0, e = a.length; e > d; d++)
                if (a[d].lat.toFixed(c) != b[d].lat.toFixed(c) || a[d].lng.toFixed(c) != b[d].lng.toFixed(c)) return !1;
            return !0
        },
        _showRoute_RedrawHistory: function(a, b) {
            if (a) {
                var c = OSRM.RoutingGeometry._decode(a.route_geometry, OSRM.C.PRECISION);
                this._history[b].route.setPositions(c), this._updateHints(a, b)
            }
        },
        _getRoute_RedrawHistory: function() {
            for (var a = 0, b = this._history_length; b > a; a++) this._history[a].markers.length > 0 && (OSRM.JSONP.clear("history" + a), OSRM.JSONP.call(this._buildCall(a) + "&instructions=false", this._callback_redrawHistory, OSRM.JSONP.empty, OSRM.DEFAULTS.JSONP_TIMEOUT, "history" + a, a))
        },
        _buildCall: function(a) {
            var b = OSRM.G.active_routing_server_url;
            b += "?z=" + OSRM.G.map.getZoom() + "&output=json&jsonp=%jsonp", this._history[a].checksum && (b += "&checksum=" + this._history[a].checksum);
            for (var c = this._history[a].markers, d = OSRM.C.PRECISION, e = 0, f = c.length; f > e; e++) b += "&loc=" + c[e].lat.toFixed(d) + "," + c[e].lng.toFixed(d), c[e].hint && (b += "&hint=" + c[e].hint);
            return b
        },
        _updateHints: function(a, b) {
            this._history[b].checksum = a.hint_data.checksum;
            for (var c = a.hint_data.locations, d = 0; d < c.length; d++) this._history[b].markers[d].hint = c[d]
        }
    }), OSRM.GUI = {
        init_functions: [],
        init: function() {
            for (var a = 0, b = OSRM.GUI.init_functions.length; b > a; a++) OSRM.GUI.init_functions[a]()
        },
        extend: function(a) {
            for (property in a) "init" == property ? OSRM.GUI.init_functions.push(a[property]) : OSRM.GUI[property] = a[property]
        }
    }, OSRM.GUIBoxGroup = function() {
        this._handles = []
    }, OSRM.extend(OSRM.GUIBoxGroup, {
        add: function(a) {
            this._handles.push(a), a.$addToGroup(this)
        },
        select: function(a) {
            for (var b = 0; b < this._handles.length; b++) this._handles[b] != a ? this._handles[b].$hideBox() : this._handles[b].$showBox()
        },
        $hide: function() {
            for (var a = 0; a < this._handles.length; a++) this._handles[a].$hide()
        },
        $show: function() {
            for (var a = 0; a < this._handles.length; a++) this._handles[a].$show()
        }
    }), OSRM.GUIBoxHandle = function(a, b, c, d, e) {
        var f = document.getElementById(a + "-toggle");
        if (null == f) return void console.log("[error] No toggle button for " + a);
        var g = document.createElement("div");
        g.id = a + "-handle-wrapper", g.className = "box-wrapper box-handle-wrapper-" + b, g.style.cssText += c;
        var h = document.createElement("div");
        h.id = a + "-handle-content", h.className = "box-content box-handle-content-" + b;
        var i = document.createElement("div");
        i.id = a + "-handle-icon", i.className = "iconic-button", i.title = a, h.appendChild(i), g.appendChild(h), document.body.appendChild(g), this._box = document.getElementById(a + "-wrapper"), this._class = this._box.className, this._width = this._box.clientWidth, this._side = b, this._handle = g, this._box_group = null, this._transitionEndFct = e, this._box.style[this._side] = -this._width + "px", this._box_visible = !1, this._box.style.visibility = "hidden", this._handle.style.visibility = "visible";
        var j = d ? OSRM.concat(this._toggle, d) : this._toggle,
            k = OSRM.bind(this, j);
        f.onclick = k, i.onclick = k;
        var j = e ? OSRM.concat(this._onTransitionEnd, e) : this._onTransitionEnd,
            k = OSRM.bind(this, j);
        if (-1 == OSRM.Browser.FF3 && -1 == OSRM.Browser.IE6_9) {
            var l = document.getElementById(a + "-wrapper");
            l.addEventListener("transitionend", k, !1), l.addEventListener("webkitTransitionEnd", k, !1), l.addEventListener("oTransitionEnd", k, !1), l.addEventListener("MSTransitionEnd", k, !1)
        } else this._legacyTransitionEndFct = OSRM.bind(this, function() {
            k({
                target: this._box
            })
        })
    }, OSRM.extend(OSRM.GUIBoxHandle, {
        boxVisible: function() {
            return this._box_visible
        },
        boxWidth: function() {
            return this._width
        },
        $addToGroup: function(a) {
            this._box_group = a
        },
        $show: function() {
            this._handle.style.visibility = "visible"
        },
        $hide: function() {
            this._handle.style.visibility = "hidden"
        },
        $showBox: function() {
            this._box_visible = !0, this._box.style.visibility = "visible", this._handle.style.visibility = "hidden", this._box.style[this._side] = "5px"
        },
        $hideBox: function() {
            this._box_visible = !1, this._box.style.visibility = "hidden", this._handle.style.visibility = "visible", this._box.style[this._side] = -this._width + "px"
        },
        _toggle: function() {
            this._box.className += " box-animated", 0 == this._box_visible ? (this._box_group.$hide(), this._box.style[this._side] = "5px", this._box.style.visibility = "visible") : this._box.style[this._side] = -this._width + "px", (-1 != OSRM.Browser.FF3 || -1 != OSRM.Browser.IE6_9) && setTimeout(this._legacyTransitionEndFct, 0)
        },
        _onTransitionEnd: function(a) {
            a.target == this._box && (this._box.className = this._class, 1 == this._box_visible ? (this._box_group.$show(), this._box_visible = !1, this._box.style.visibility = "hidden") : (this._box_visible = !0, this._box.style.visibility = "visible"))
        }
    }), OSRM.GUI.extend({
        selectorInit: function(a, b, c, d) {
            var e = document.getElementById(a);
            e.className += " styled-select-helper base-font", e.onchange = function() {
                OSRM.GUI._selectorOnChange(this), d(this.value)
            };
            for (var f = 0, g = b.length; g > f; f++) {
                var h = document.createElement("option");
                h.innerHTML = b[f].display, h.value = b[f].value, e.appendChild(h)
            }
            e.value = b[c].value;
            var i = document.createTextNode(b[c].display),
                j = document.createElement("span");
            j.className = "styled-select base-font", j.id = "styled-select-" + e.id, j.appendChild(i), e.parentNode.insertBefore(j, e), j.style.width = e.offsetWidth - 2 + "px", j.style.height = e.offsetHeight + "px"
        },
        _selectorOnChange: function(a) {
            for (var b = a.getElementsByTagName("option"), c = 0; c < b.length; c++)
                if (1 == b[c].selected) {
                    document.getElementById("styled-select-" + a.id).childNodes[0].nodeValue = b[c].childNodes[0].nodeValue;
                    break
                }
        },
        selectorChange: function(a, b) {
            var c = document.getElementById(a);
            c.value = b, OSRM.GUI._selectorOnChange(c)
        },
        selectorRenameOptions: function(a, b) {
            var c = document.getElementById(a),
                d = document.getElementById("styled-select-" + a),
                e = document.createElement("select");
            e.id = a, e.className = c.className, e.onchange = c.onchange;
            for (var f = "", g = 0, h = b.length; h > g; g++) {
                var i = document.createElement("option");
                i.innerHTML = b[g].display, i.value = b[g].value, e.appendChild(i), b[g].value == c.value && (f = b[g].display)
            }
            e.value = c.value, c.parentNode.insertBefore(e, c), c.parentNode.removeChild(c), d.childNodes[0].nodeValue = f, d.style.width = e.offsetWidth - 2 + "px", d.style.height = e.offsetHeight + "px"
        }
    }), OSRM.GUI.extend({
        init: function() {
            var a = new OSRM.GUIBoxGroup;
            OSRM.G.main_handle = new OSRM.GUIBoxHandle("main", "left", "left:-5px;top:5px;", OSRM.GUI.beforeMainTransition, OSRM.GUI.afterMainTransition), a.add(OSRM.G.main_handle), a.select(OSRM.G.main_handle);
            var b = new OSRM.GUIBoxGroup,
                c = new OSRM.GUIBoxHandle("config", "right", "right:-5px;bottom:70px;"),
                d = new OSRM.GUIBoxHandle("mapping", "right", "right:-5px;bottom:25px;");
            b.add(c), b.add(d), b.select(null), document.getElementById("gui-input-source").value = OSRM.DEFAULTS.ONLOAD_SOURCE, document.getElementById("gui-input-target").value = OSRM.DEFAULTS.ONLOAD_TARGET, OSRM.GUI.initDistanceFormatsSelector()
        },
        setLabels: function() {
            document.getElementById("open-editor").innerHTML = OSRM.loc("OPEN_EDITOR"), document.getElementById("open-josm").innerHTML = OSRM.loc("OPEN_JOSM"), document.getElementById("open-osmbugs").innerHTML = OSRM.loc("OPEN_OSMBUGS"), document.getElementById("gui-reset").innerHTML = OSRM.loc("GUI_RESET"), document.getElementById("gui-reverse").innerHTML = OSRM.loc("GUI_REVERSE"), document.getElementById("gui-option-highlight-nonames-label").lastChild.nodeValue = OSRM.loc("GUI_HIGHLIGHT_UNNAMED_ROADS"), document.getElementById("gui-option-show-previous-routes-label").lastChild.nodeValue = OSRM.loc("GUI_SHOW_PREVIOUS_ROUTES"), document.getElementById("gui-search-source").innerHTML = OSRM.loc("GUI_SEARCH"), document.getElementById("gui-search-target").innerHTML = OSRM.loc("GUI_SEARCH"), document.getElementById("gui-search-source-label").innerHTML = OSRM.loc("GUI_START") + ":", document.getElementById("gui-search-target-label").innerHTML = OSRM.loc("GUI_END") + ":", document.getElementById("gui-input-source").title = OSRM.loc("GUI_START_TOOLTIP"), document.getElementById("gui-input-target").title = OSRM.loc("GUI_END_TOOLTIP"), document.getElementById("legal-notice").innerHTML = OSRM.loc("GUI_LEGAL_NOTICE"), document.getElementById("gui-mapping-label").innerHTML = OSRM.loc("GUI_MAPPING_TOOLS"), document.getElementById("gui-external-tools-label").innerHTML = OSRM.loc("GUI_EXTERNAL_TOOLS"), document.getElementById("gui-config-label").innerHTML = OSRM.loc("GUI_CONFIGURATION"), document.getElementById("gui-language-2-label").innerHTML = OSRM.loc("GUI_LANGUAGE") + ":", document.getElementById("gui-units-label").innerHTML = OSRM.loc("GUI_UNITS") + ":", document.getElementById("gui-data-timestamp-label").innerHTML = OSRM.loc("GUI_DATA_TIMESTAMP"), document.getElementById("gui-data-timestamp").innerHTML = OSRM.G.data_timestamp, document.getElementById("gui-timestamp-label").innerHTML = OSRM.loc("GUI_VERSION"), document.getElementById("gui-timestamp").innerHTML = OSRM.DATE + "; v" + OSRM.VERSION, document.getElementById("config-handle-icon").title = OSRM.loc("GUI_CONFIGURATION"), document.getElementById("mapping-handle-icon").title = OSRM.loc("GUI_MAPPING_TOOLS"), document.getElementById("main-handle-icon").title = OSRM.loc("GUI_MAIN_WINDOW"), OSRM.G.map.zoomControl.setTooltips(OSRM.loc("GUI_ZOOM_IN"), OSRM.loc("GUI_ZOOM_OUT")), OSRM.G.map.locationsControl.setTooltips(OSRM.loc("GUI_ZOOM_ON_USER"), OSRM.loc("GUI_ZOOM_ON_ROUTE")), OSRM.GUI.setDistanceFormatsLanguage(), OSRM.GUI.setRoutingEnginesLanguage()
        },
        clearResults: function() {
            document.getElementById("information-box").className = "information-box-with-normal-header", document.getElementById("information-box").innerHTML = "", document.getElementById("information-box-header").innerHTML = ""
        },
        beforeMainTransition: function() {
            OSRM.G.map.zoomControl.hide()
        },
        afterMainTransition: function() {
            OSRM.G.map.zoomControl.show()
        },
        initDistanceFormatsSelector: function() {
            var a = OSRM.GUI.getDistanceFormats();
            OSRM.GUI.selectorInit("gui-units-toggle", a, OSRM.DEFAULTS.DISTANCE_FORMAT, OSRM.GUI._onDistanceFormatChanged)
        },
        setDistanceFormat: function(a) {
            OSRM.G.active_distance_format != a && (OSRM.G.active_distance_format = a, OSRM.G.map && (OSRM.G.map.scaleControl.removeFrom(OSRM.G.map), OSRM.G.map.scaleControl.options.metric = 1 != a, OSRM.G.map.scaleControl.options.imperial = 1 == a, OSRM.G.map.scaleControl.addTo(OSRM.G.map)), OSRM.Utils.toHumanDistance = 1 == a ? OSRM.Utils.toHumanDistanceMiles : OSRM.Utils.toHumanDistanceMeters)
        },
        _onDistanceFormatChanged: function(a) {
            OSRM.GUI.setDistanceFormat(a), OSRM.Routing.getRoute({
                keepAlternative: !0
            })
        },
        setDistanceFormatsLanguage: function() {
            var a = OSRM.GUI.getDistanceFormats();
            OSRM.GUI.selectorRenameOptions("gui-units-toggle", a)
        },
        getDistanceFormats: function() {
            return [{
                display: OSRM.loc("GUI_KILOMETERS"),
                value: 0
            }, {
                display: OSRM.loc("GUI_MILES"),
                value: 1
            }]
        },
        queryDataTimestamp: function() {
            OSRM.G.data_timestamp = "n/a", document.getElementById("gui-data-timestamp").innerHTML = OSRM.G.data_timestamp, OSRM.JSONP.call(OSRM.G.active_routing_timestamp_url + "?jsonp=%jsonp", OSRM.GUI.setDataTimestamp, OSRM.JSONP.empty, OSRM.DEFAULTS.JSONP_TIMEOUT, "data_timestamp")
        },
        setDataTimestamp: function(a) {
            a && (OSRM.G.data_timestamp = a.timestamp.slice(0, 25).replace(/<\/?[^>]+(>|$)/g, ""), document.getElementById("gui-data-timestamp").innerHTML = OSRM.G.data_timestamp)
        }
    }), OSRM.GUI.extend({
        activeExclusive: void 0,
        activeTooltip: void 0,
        tooltips: {},
        init: function() {
            var a = OSRM.DEFAULTS.NOTIFICATIONS;
            if (1 == a.MAINTENANCE) {
                var b = OSRM.DEFAULTS.OVERRIDE_MAINTENANCE_NOTIFICATION_HEADER || OSRM.loc("NOTIFICATION_MAINTENANCE_HEADER"),
                    c = OSRM.DEFAULTS.OVERRIDE_MAINTENANCE_NOTIFICATION_BODY || OSRM.loc("NOTIFICATION_MAINTENANCE_BODY");
                return OSRM.GUI.exclusiveNotify(b, c, !1), void(OSRM.GUI.activeExclusive = "MAINTENANCE")
            }
            var a = OSRM.DEFAULTS.NOTIFICATIONS,
                d = OSRM.GUI.tooltips;
            for (id in a) OSRM.Utils.isNumber(a[id]) && (d[id] = {}, d[id]._timer = setTimeout(function(a) {
                return function() {
                    OSRM.GUI._showTooltip(a)
                }
            }(id), a[id]), d[id]._pending = !0)
        },
        deactivateTooltip: function(a) {
            var b = OSRM.GUI.tooltips;
            void 0 != b[a] && (b[a]._pending = !1)
        },
        deactivateTooltips: function() {
            var a = OSRM.GUI.tooltips;
            for (id in a) a[id]._pending = !1
        },
        _showTooltip: function(a) {
            var b = OSRM.GUI.tooltips;
            if (void 0 != b[a] && 0 != b[a]._pending) {
                var c = OSRM.DEFAULTS.NOTIFICATIONS;
                if (OSRM.GUI.isTooltipVisible()) return void(b[a]._timer = setTimeout(function(a) {
                    return function() {
                        OSRM.GUI._showTooltip(a)
                    }
                }(a), c[a]));
                OSRM.GUI.tooltipNotify(OSRM.loc("NOTIFICATION_" + a + "_HEADER"), OSRM.loc("NOTIFICATION_" + a + "_BODY")), OSRM.GUI.activeTooltip = a, b[a]._pending = !1
            }
        },
        exclusiveNotify: function(a, b, c) {
            document.getElementById("exclusive-notification-blanket").style.display = "block", document.getElementById("exclusive-notification-label").innerHTML = a, document.getElementById("exclusive-notification-box").innerHTML = b, c ? (document.getElementById("exclusive-notification-toggle").style.display = "block", document.getElementById("exclusive-notification-toggle").onclick = OSRM.GUI.exclusiveDenotify) : document.getElementById("exclusive-notification-toggle").style.display = "none", OSRM.GUI.exclusiveResize()
        },
        exclusiveDenotify: function() {
            document.getElementById("exclusive-notification-blanket").style.display = "none", OSRM.GUI.activeExclusive = void 0
        },
        exclusiveUpdate: function() {
            if (void 0 != OSRM.GUI.activeExclusive) {
                var a = OSRM.DEFAULTS["OVERRIDE_" + OSRM.GUI.activeExclusive + "_HEADER"] || OSRM.loc("NOTIFICATION_MAINTENANCE_HEADER"),
                    b = OSRM.DEFAULTS["OVERRIDE_" + OSRM.GUI.activeExclusive + "_BODY"] || OSRM.loc("NOTIFICATION_MAINTENANCE_BODY");
                document.getElementById("exclusive-notification-label").innerHTML = a, document.getElementById("exclusive-notification-box").innerHTML = b, OSRM.GUI.exclusiveResize()
            }
        },
        exclusiveResize: function() {
            var a = document.getElementById("exclusive-notification-box").clientHeight;
            document.getElementById("exclusive-notification-content").style.height = a + 28 + "px", document.getElementById("exclusive-notification-wrapper").style.height = a + 48 + "px"
        },
        inMaintenance: function() {
            return "MAINTENANCE" == OSRM.GUI.activeExclusive
        },
        tooltipNotify: function(a, b) {
            document.getElementById("tooltip-notification-wrapper").style.display = "block", document.getElementById("tooltip-notification-label").innerHTML = a, document.getElementById("tooltip-notification-box").innerHTML = b, document.getElementById("tooltip-notification-box").style.display = "block", OSRM.GUI.tooltipResize(), document.getElementById("tooltip-notification-toggle").onclick = OSRM.GUI.tooltipDenotify, document.getElementById("tooltip-notification-resize").onclick = OSRM.GUI.tooltipResize
        },
        tooltipDenotify: function() {
            document.getElementById("tooltip-notification-wrapper").style.display = "none", OSRM.GUI.activeTooltip = void 0
        },
        tooltipUpdate: function() {
            void 0 != OSRM.GUI.activeTooltip && (document.getElementById("tooltip-notification-label").innerHTML = OSRM.loc("NOTIFICATION_" + OSRM.GUI.activeTooltip + "_HEADER"), document.getElementById("tooltip-notification-box").innerHTML = OSRM.loc("NOTIFICATION_" + OSRM.GUI.activeTooltip + "_BODY"), OSRM.GUI.tooltipResize(), OSRM.GUI.tooltipResize())
        },
        tooltipResize: function() {
            if ("none" == document.getElementById("tooltip-notification-box").style.display) {
                document.getElementById("tooltip-notification-box").style.display = "block";
                var a = document.getElementById("tooltip-notification-box").clientHeight;
                document.getElementById("tooltip-notification-content").style.height = a + 28 + "px", document.getElementById("tooltip-notification-wrapper").style.height = a + 48 + "px", document.getElementById("tooltip-notification-resize").className = "iconic-button up-marker top-right-button"
            } else document.getElementById("tooltip-notification-box").style.display = "none", document.getElementById("tooltip-notification-content").style.height = "18px", document.getElementById("tooltip-notification-wrapper").style.height = "38px", document.getElementById("tooltip-notification-resize").className = "iconic-button down-marker top-right-button"
        },
        isTooltipVisible: function() {
            return "block" == document.getElementById("tooltip-notification-wrapper").style.display
        },
        updateNotifications: function() {
            OSRM.GUI.exclusiveUpdate(), OSRM.GUI.tooltipUpdate()
        }
    }), OSRM.GLOBALS.route = null, OSRM.GLOBALS.markers = null, OSRM.GLOBALS.dragging = null, OSRM.GLOBALS.dragid = null, OSRM.GLOBALS.pending = !1, OSRM.Routing = {
        init: function() {
            OSRM.GUI.setRoutingEngine(OSRM.DEFAULTS.ROUTING_ENGINE), OSRM.G.markers = new OSRM.Markers, OSRM.G.route = new OSRM.Route, OSRM.G.response = {
                via_points: []
            }, OSRM.RoutingDescription.init()
        },
        timeoutRoute: function() {
            OSRM.G.response = {
                via_points: []
            }, OSRM.RoutingGeometry.showNA(), OSRM.RoutingNoNames.showNA(), OSRM.RoutingDescription.showNA(OSRM.loc("TIMED_OUT")), OSRM.Routing._snapRoute()
        },
        timeoutRoute_Dragging: function() {
            OSRM.RoutingGeometry.showNA(), OSRM.RoutingDescription.showNA(OSRM.loc("TIMED_OUT"))
        },
        timeoutRoute_Reversed: function() {
            OSRM.G.markers.reverseMarkers(), OSRM.Routing.timeoutRoute()
        },
        showRoute: function(a, b) {
            if (a && (1 != b.keepAlternative && (OSRM.G.active_alternative = 0), OSRM.G.response = a, OSRM.Routing._snapRoute(), 207 == a.status ? (OSRM.RoutingGeometry.showNA(), OSRM.RoutingNoNames.showNA(), OSRM.RoutingDescription.showNA(OSRM.loc("NO_ROUTE_FOUND"))) : (("undefined" == typeof a.found_alternative || a.found_alternative) && OSRM.RoutingAlternatives.prepare(OSRM.G.response), OSRM.RoutingGeometry.show(OSRM.G.response), OSRM.RoutingNoNames.show(OSRM.G.response), OSRM.RoutingDescription.show(OSRM.G.response)), OSRM.Routing._updateHints(a), 1 == b.recenter)) {
                var c = new L.LatLngBounds(OSRM.G.route._current_route.getPositions());
                OSRM.G.map.setViewBoundsUI(c)
            }
        },
        showRoute_Dragging: function(a) {
            a && OSRM.G.dragging && (OSRM.G.response = a, 207 == a.status ? (OSRM.RoutingGeometry.showNA(), OSRM.RoutingDescription.showNA(OSRM.loc("YOUR_ROUTE_IS_BEING_COMPUTED"))) : (OSRM.RoutingGeometry.show(a), OSRM.RoutingDescription.showSimple(a)), OSRM.Routing._updateHints(a), OSRM.G.pending && setTimeout(OSRM.Routing.draggingTimeout, 1))
        },
        showRoute_Redraw: function(a, b) {
            a && (0 == b.keepAlternative && (OSRM.G.active_alternative = 0), OSRM.G.response = a, 207 != a.status && (("undefined" == typeof a.found_alternative || a.found_alternative) && OSRM.RoutingAlternatives.prepare(OSRM.G.response), OSRM.RoutingGeometry.show(OSRM.G.response), OSRM.RoutingNoNames.show(OSRM.G.response)), OSRM.Routing._updateHints(a))
        },
		timeoutTables: function(data){
			document.getElementById('time-table-box').innerHTML='';
			document.getElementById('length-table-box').innerHTML='';
		},
		showTables: function(data)
		{
			var n=data.time_table[0].length;
			var t='<b>Матрица времени</b><table><tr><th>&nbsp;</th>';
			for(var i=0;i<n;++i) t+='<th>'+(i==0?'start':i==n-1?'end':i)+'</th>';
			t+='</tr>';
			for(var i=0;i<n;++i) {
				t+='<tr><th>'+(i==0?'start':i==n-1?'end':i)+'</th>';
				for(var j=0;j<n;++j) 
					t+='<td>'+data.time_table[0][i][j]/10+'</td>';
				t+='</tr>';
			}
			t+='</table>';
			document.getElementById('time-table-box').innerHTML=t;
			var t='<b>Матрица расстояния</b><table><tr><th>&nbsp;</th>';
			for(var i=0;i<n;++i) t+='<th>'+(i==0?'start':i==n-1?'end':i)+'</th>';
			t+='</tr>';
			for(var i=0;i<n;++i) {
				t+='<tr><th>'+(i==0?'start':i==n-1?'end':i)+'</th>';
				for(var j=0;j<n;++j) 
					t+='<td>'+data.length_table[0][i][j]+'</td>';
				t+='</tr>';
			}
			t+='</table>';
			document.getElementById('length-table-box').innerHTML=t;
		},
        getRoute: function(a) {
			//OSRM.JSONP.call(OSRM.Routing._buildCall()+'&compression=false', function(a){console.log(a.route_geometry)}, function(){}, OSRM.DEFAULTS.JSONP_TIMEOUT, "log");
			if(OSRM.G.markers.route.length > 1) OSRM.JSONP.call(OSRM.Routing._buildCall().replace('viaroute','math'), showChains, function(){}, OSRM.DEFAULTS.JSONP_TIMEOUT*100, "math")
            return OSRM.G.markers.route.length < 2 ? (void OSRM.G.route.hideRoute(), OSRM.Routing.timeoutTables()) : (a = a || {}, OSRM.JSONP.clear("dragging"), OSRM.JSONP.clear("redraw"), OSRM.JSONP.clear("table"), OSRM.JSONP.clear("route"), void OSRM.JSONP.call(OSRM.Routing._buildCall() + "&instructions=true", OSRM.Routing.showRoute, OSRM.Routing.timeoutRoute, OSRM.DEFAULTS.JSONP_TIMEOUT, "route", a), OSRM.JSONP.call(OSRM.Routing._buildCall().replace('viaroute','table'), OSRM.Routing.showTables, OSRM.Routing.timeoutTables, OSRM.DEFAULTS.JSONP_TIMEOUT, "table"))
        },
        getRoute_Reversed: function() {
            OSRM.G.markers.route.length < 2 || (OSRM.JSONP.clear("dragging"), OSRM.JSONP.clear("redraw"), OSRM.JSONP.clear("route"), OSRM.JSONP.call(OSRM.Routing._buildCall() + "&instructions=true", OSRM.Routing.showRoute, OSRM.Routing.timeoutRoute_Reversed, OSRM.DEFAULTS.JSONP_TIMEOUT, "route", {}))
        },
        getRoute_Redraw: function(a) {
            OSRM.G.markers.route.length < 2 || (a = a || {}, OSRM.JSONP.clear("dragging"), OSRM.JSONP.clear("redraw"), OSRM.JSONP.call(OSRM.Routing._buildCall() + "&instructions=true", OSRM.Routing.showRoute_Redraw, OSRM.Routing.timeoutRoute, OSRM.DEFAULTS.JSONP_TIMEOUT, "redraw", a))
        },
        getRoute_Dragging: function() {
            OSRM.G.pending = !OSRM.JSONP.call(OSRM.Routing._buildCall() + "&alt=false&instructions=false", OSRM.Routing.showRoute_Dragging, OSRM.Routing.timeoutRoute_Dragging, OSRM.DEFAULTS.JSONP_TIMEOUT, "dragging"), OSRM.JSONP.call(OSRM.Routing._buildCall().replace('viaroute','table'), OSRM.Routing.showTables, OSRM.Routing.timeoutTables, OSRM.DEFAULTS.JSONP_TIMEOUT, "table")
        },
        draggingTimeout: function() {
            OSRM.G.markers.route[OSRM.G.dragid].hint = null, OSRM.Routing.getRoute_Dragging()
        },
        _buildCall: function() {
            var a = OSRM.G.active_routing_server_url;
            a += "?z=" + OSRM.G.map.getZoom() + "&output=json&jsonp=%jsonp", OSRM.G.markers.checksum && (a += "&checksum=" + OSRM.G.markers.checksum);
            for (var b = OSRM.G.markers.route, c = OSRM.C.PRECISION, d = 0, e = b.length; e > d; d++) a += "&loc=" + b[d].getLat().toFixed(c) + "," + b[d].getLng().toFixed(c), b[d].hint && (a += "&hint=" + b[d].hint);
			var load=document.getElementById('gui-input-load').value;
			if(!load) load=0;
			var height=document.getElementById('gui-input-height').value;
			if(!height) height=0;
			a+='&transport='+load+','+height;
            return a
        },
        _updateHints: function(a) {
            var b = a.hint_data.locations;
            OSRM.G.markers.checksum = a.hint_data.checksum;
            for (var c = 0; c < b.length; c++) OSRM.G.markers.route[c].hint = b[c]
        },
        _snapRoute: function() {
            for (var a = OSRM.G.markers.route, b = OSRM.G.response.via_points, c = 0; c < b.length; c++) a[c].setPosition(new L.LatLng(b[c][0], b[c][1]));
            OSRM.Geocoder.updateAddress(OSRM.C.SOURCE_LABEL), OSRM.Geocoder.updateAddress(OSRM.C.TARGET_LABEL), OSRM.G.markers.relabelViaMarkers()
        }
    }, OSRM.RoutingAlternatives = {
        _buttons: [{
            id: "gui-a",
            label: "A"
        }, {
            id: "gui-b",
            label: "B"
        }],
        init: function() {
            OSRM.G.active_alternative = 0, OSRM.G.alternative_count = 0
        },
        prepare: function(a) {
            var b = OSRM.G.response;
            b.route_name = b.route_name || [], b.alternative_names = b.alternative_names || [
                []
            ], b.alternative_geometries && (b.alternative_geometries.unshift(a.route_geometry), b.alternative_instructions.unshift(a.route_instructions), b.alternative_summaries.unshift(a.route_summary), b.alternative_names.unshift(a.route_name), OSRM.G.alternative_count = a.alternative_geometries.length, OSRM.G.active_alternative >= OSRM.G.alternative_count && (OSRM.G.active_alternative = 0), b.route_geometry = b.alternative_geometries[OSRM.G.active_alternative], b.route_instructions = b.alternative_instructions[OSRM.G.active_alternative], b.route_summary = b.alternative_summaries[OSRM.G.active_alternative], b.route_name = b.alternative_names[OSRM.G.active_alternative])
        },
        setActive: function(a) {
            OSRM.G.active_alternative = a;
            for (var b = OSRM.RoutingAlternatives._buttons, c = 0, d = OSRM.G.alternative_count; d > c; c++) document.getElementById(b[c].id).className = a == c ? "button-pressed top-right-button" : "button top-right-button"
        },
        show: function() {
            for (var a = OSRM.RoutingAlternatives._buttons, b = "", c = 0, d = OSRM.G.alternative_count; d > c && OSRM.G.response.alternative_summaries; c++) {
                for (var e = OSRM.Utils.toHumanDistance(OSRM.G.response.alternative_summaries[c].total_distance), f = OSRM.Utils.toHumanTime(OSRM.G.response.alternative_summaries[c].total_time), g = " &#10;(", h = 0, i = OSRM.G.response.alternative_names[c].length; i > h; h++) g += (h > 0 && "" != OSRM.G.response.alternative_names[c][h] && "" != OSRM.G.response.alternative_names[c][h - 1] ? " - " : "") + OSRM.G.response.alternative_names[c][h];
                g += ")";
                var j = OSRM.loc("DISTANCE") + ": " + e + " &#10;" + OSRM.loc("DURATION") + ": " + f + g,
                    k = c == OSRM.G.active_alternative ? "button-pressed" : "button";
                b = '<a class="' + k + ' top-right-button" id="' + a[c].id + '" title="' + j + '">' + a[c].label + "</a>" + b
            }
            for (var c = OSRM.G.alternative_count, d = a.length; d > c; ++c) b = '<a class="button-inactive top-right-button" id="' + a[c].id + '">' + a[c].label + "</a>" + b;
            document.getElementById("information-box-header").innerHTML = b + document.getElementById("information-box-header").innerHTML;
            for (var c = 0, d = OSRM.G.alternative_count; d > c; c++) document.getElementById(a[c].id).onclick = function(a) {
                return function() {
                    OSRM.RoutingAlternatives._click(a)
                }
            }(c), document.getElementById(a[c].id).onmouseover = function(a) {
                return function() {
                    OSRM.RoutingAlternatives._mouseover(a)
                }
            }(c), document.getElementById(a[c].id).onmouseout = function(a) {
                return function() {
                    OSRM.RoutingAlternatives._mouseout(a)
                }
            }(c)
        },
        _click: function(a) {
            if (OSRM.G.active_alternative != a) {
                OSRM.RoutingAlternatives.setActive(a), OSRM.G.route.hideAlternativeRoute();
                var b = OSRM.G.response;
                b.route_geometry = b.alternative_geometries[a], b.route_instructions = b.alternative_instructions[a], b.route_summary = b.alternative_summaries[a], b.route_name = b.alternative_names[a], OSRM.RoutingGeometry.show(b), OSRM.RoutingNoNames.show(b), OSRM.RoutingDescription.show(b), OSRM.G.markers.highlight.hide()
            }
        },
        _mouseover: function(a) {
            if (OSRM.G.active_alternative != a) {
                var b = OSRM.RoutingGeometry._decode(OSRM.G.response.alternative_geometries[a], OSRM.C.PRECISION);
                OSRM.G.route.showAlternativeRoute(b)
            }
        },
        _mouseout: function(a) {
            OSRM.G.active_alternative != a && OSRM.G.route.hideAlternativeRoute()
        }
    }, OSRM.RoutingDescription = {
        QR_DIRECTORY: "qrcodes/",
        init: function() {
            OSRM.G.active_shortlink = null, OSRM.Browser.onUnloadHandler(OSRM.RoutingDescription.uninit)
        },
        uninit: function() {
            OSRM.G.qrcodewindow && OSRM.G.qrcodewindow.close()
        },
        onMouseOverRouteDescription: function(a, b) {
            OSRM.G.markers.hover.setPosition(new L.LatLng(a, b)), OSRM.G.markers.hover.show()
        },
        onMouseOutRouteDescription: function() {
            OSRM.G.markers.hover.hide()
        },
        onClickRouteDescription: function(a, b, c) {
            OSRM.G.markers.highlight.setPosition(new L.LatLng(a, b)), OSRM.G.markers.highlight.show(), OSRM.G.markers.highlight.centerView(OSRM.DEFAULTS.HIGHLIGHT_ZOOM_LEVEL), null != OSRM.G.markers.highlight.description && document.getElementById("description-" + OSRM.G.markers.highlight.description) && (document.getElementById("description-" + OSRM.G.markers.highlight.description).className = "description-body-item"), OSRM.G.markers.highlight.description = c, document.getElementById("description-" + c).className = "description-body-item description-body-item-selected"
        },
        onClickCreateShortcut: function(a) {
            var b = OSRM.C.PRECISION;
            if (a += "&z=" + OSRM.G.map.getZoom() + "&center=" + OSRM.G.map.getCenter().lat.toFixed(b) + "," + OSRM.G.map.getCenter().lng.toFixed(b), a += "&alt=" + OSRM.G.active_alternative, a += "&df=" + OSRM.G.active_distance_format, a += "&re=" + OSRM.G.active_routing_engine, a += "&ly=" + OSRM.Utils.getHash(OSRM.G.map.layerControl.getActiveLayerName()), "" == OSRM.DEFAULTS.HOST_SHORTENER_URL) {
                var c = {};
                return c.label = OSRM.loc("ROUTE_LINK"), c[OSRM.DEFAULTS.SHORTENER_REPLY_PARAMETER] = a, void OSRM.RoutingDescription.showRouteLink(c)
            }
            var d = OSRM.DEFAULTS.HOST_SHORTENER_URL + OSRM.DEFAULTS.SHORTENER_PARAMETERS.replace(/%url/, a);
            OSRM.JSONP.call(d, OSRM.RoutingDescription.showRouteLink, OSRM.RoutingDescription.showRouteLink_TimeOut, OSRM.DEFAULTS.JSONP_TIMEOUT, "shortener"), document.getElementById("route-link").innerHTML = '[<a class="text-link-inactive">' + OSRM.loc("GENERATE_LINK_TO_ROUTE") + "</a>]"
        },
        showRouteLink: function(a) {
            if (!a || !a[OSRM.DEFAULTS.SHORTENER_REPLY_PARAMETER]) return void OSRM.RoutingDescription.showRouteLink_TimeOut();
            OSRM.G.active_shortlink = a[OSRM.DEFAULTS.SHORTENER_REPLY_PARAMETER];
            var b = a.label;
            b || (b = OSRM.G.active_shortlink.substring(7)), document.getElementById("route-link").innerHTML = '[<a class="text-link" onClick="OSRM.RoutingDescription.showQRCode();">' + OSRM.loc("QR") + '</a>]&nbsp;[<a class="text-link" href="' + OSRM.G.active_shortlink + '">' + b + "</a>]"
        },
        showRouteLink_TimeOut: function() {
            document.getElementById("route-link").innerHTML = '[<a class="text-link-inactive">' + OSRM.loc("LINK_TO_ROUTE_TIMEOUT") + "</a>]"
        },
        showQRCode: function() {
            OSRM.G.qrcodewindow && OSRM.G.qrcodewindow.close(), OSRM.G.qrcodewindow = window.open(OSRM.RoutingDescription.QR_DIRECTORY + "qrcodes.html", "", "width=280,height=250,left=100,top=100,dependent=yes,location=no,menubar=no,scrollbars=no,status=no,toolbar=no,resizable=no")
        },
        show: function(a) {
            var b = OSRM.C.PRECISION;
            OSRM.GUI.activateRouteFeatures();
            for (var c = "?hl=" + OSRM.Localization.current_language, d = 0; d < OSRM.G.markers.route.length; d++) c += "&loc=" + OSRM.G.markers.route[d].getLat().toFixed(b) + "," + OSRM.G.markers.route[d].getLng().toFixed(b);
            var e = '[<a class="text-link" onclick="OSRM.RoutingDescription.onClickCreateShortcut(\'' + OSRM.DEFAULTS.WEBSITE_URL + c + "')\">" + OSRM.loc("GET_LINK_TO_ROUTE") + "</a>]",
                f = '[<a class="text-link" onClick="document.location.href=\'' + OSRM.G.active_routing_server_url + c + "&output=gpx';\">" + OSRM.loc("GPX_FILE") + "</a>]",
                g = null;
            OSRM.G.markers.highlight.isShown() && (g = OSRM.G.markers.highlight.description);
            var h = OSRM.G.route.getPositions(),
                i = "";
            i += '<table class="description medium-font">';
            for (var d = 0; d < a.route_instructions.length; d++) {
                var j = "description-body-odd";
                d % 2 == 0 && (j = "description-body-even"), i += '<tr class="' + j + '">', i += '<td class="description-body-directions">', i += '<img class="description-body-direction" src="' + OSRM.RoutingDescription._getDrivingInstructionIcon(a.route_instructions[d][0]) + '" alt=""/>', i += "</td>", i += '<td class="description-body-items">';
                var k = h[a.route_instructions[d][3]];
                i += '<div id="description-' + d + '" class="description-body-item ' + (g == d ? "description-body-item-selected" : "") + '" onclick="OSRM.RoutingDescription.onClickRouteDescription(' + k.lat.toFixed(b) + "," + k.lng.toFixed(b) + "," + d + ')" onmouseover="OSRM.RoutingDescription.onMouseOverRouteDescription(' + k.lat.toFixed(b) + "," + k.lng.toFixed(b) + ')" onmouseout="OSRM.RoutingDescription.onMouseOutRouteDescription(' + k.lat.toFixed(b) + "," + k.lng.toFixed(b) + ')">', i += "" != a.route_instructions[d][1] ? OSRM.loc(OSRM.RoutingDescription._getDrivingInstruction(a.route_instructions[d][0])).replace(/\[(.*)\]/, "$1").replace(/%s/, OSRM.RoutingDescription._getStreetName(a.route_instructions[d][1])).replace(/%d/, OSRM.loc(a.route_instructions[d][6])) : OSRM.loc(OSRM.RoutingDescription._getDrivingInstruction(a.route_instructions[d][0])).replace(/\[(.*)\]/, "").replace(/%d/, OSRM.loc(a.route_instructions[d][6])), i += "</div>", i += "</td>", i += '<td class="description-body-distance">', d != a.route_instructions.length - 1 && (i += "<b>" + OSRM.Utils.toHumanDistance(a.route_instructions[d][2]) + "</b>"), i += "</td>", i += "</tr>"
            }
            i += "</table>";
            for (var l = "(", m = 0, n = a.route_name.length; n > m; m++) l += (m > 0 && "" != a.route_name[m] && "" != a.route_name[m - 1] ? " - " : "") + "<span style='white-space:nowrap;'>" + OSRM.RoutingDescription._getStreetName(a.route_name[m]) + "</span>";
            "(" == l && (l += " - "), l += ")", header = OSRM.RoutingDescription._buildHeader(OSRM.Utils.toHumanDistance(a.route_summary.total_distance), OSRM.Utils.toHumanTime(a.route_summary.total_time), e, f, l);
            var o = document.createElement("tempDiv");
            document.body.appendChild(o), o.className = "base-font absolute-hidden", o.innerHTML = l;
            var p = o.clientWidth,
                q = 370;
            document.body.removeChild(o), document.getElementById("information-box").className = p > q ? "information-box-with-larger-header" : "information-box-with-large-header", document.getElementById("information-box-header").innerHTML = header, document.getElementById("information-box").innerHTML = i, OSRM.RoutingAlternatives.show()
        },
        showSimple: function(a) {
            header = OSRM.RoutingDescription._buildHeader(OSRM.Utils.toHumanDistance(a.route_summary.total_distance), OSRM.Utils.toHumanTime(a.route_summary.total_time), "", ""), document.getElementById("information-box").className = "information-box-with-normal-header", document.getElementById("information-box-header").innerHTML = header, document.getElementById("information-box").innerHTML = "<div class='no-results big-font'>" + OSRM.loc("YOUR_ROUTE_IS_BEING_COMPUTED") + "</div>"
        },
        showNA: function(a) {
            OSRM.GUI.activateRouteFeatures();
            for (var b = "?hl=" + OSRM.Localization.current_language, c = OSRM.C.PRECISION, d = 0; d < OSRM.G.markers.route.length; d++) b += "&loc=" + OSRM.G.markers.route[d].getLat().toFixed(c) + "," + OSRM.G.markers.route[d].getLng().toFixed(c);
            var e = '[<a class="text-link" onclick="OSRM.RoutingDescription.onClickCreateShortcut(\'' + OSRM.DEFAULTS.WEBSITE_URL + b + "')\">" + OSRM.loc("GET_LINK_TO_ROUTE") + "</a>]";
            header = OSRM.RoutingDescription._buildHeader("N/A", "N/A", e, ""), document.getElementById("information-box").className = "information-box-with-normal-header", document.getElementById("information-box-header").innerHTML = header, document.getElementById("information-box").innerHTML = "<div class='no-results big-font'>" + a + "</div>"
        },
        _buildHeader: function(a, b, c, d, e) {
            var f = '<div class="header-title">' + OSRM.loc("ROUTE_DESCRIPTION") + (e ? '<br/><div class="header-subtitle">' + e + "</div>" : "") + '</div><div class="full"><div class="row"><div class="left"><div class="full"><div class="row"><div class="left header-label nowrap">' + OSRM.loc("DISTANCE") + ':</div><div class="left header-content stretch">' + a + '</div></div><div class="row"><div class="left header-label nowrap">' + OSRM.loc("DURATION") + ':</div><div class="left header-content stretch">' + b + '</div></div></div></div><div class="left"><div class="full"><div class="row"><div class="right header-content" id="route-link">' + c + '</div></div><div class="row"><div class="right header-content">' + d + "</div></div></div></div></div></div></div>";
            return f
        },
        _getDrivingInstructionIcon: function(a) {
            var b = "direction_";
            return a = a.replace(/^11-\d{1,}$/, "11"), b += a, OSRM.G.images[b] ? OSRM.G.images[b].getAttribute("src") : OSRM.G.images.direction_0.getAttribute("src")
        },
        _getDrivingInstruction: function(a) {
            var b = "DIRECTION_";
            a = a.replace(/^11-\d{2,}$/, "11-x"), b += a;
            var c = OSRM.loc(b);
            return c == b ? OSRM.loc("DIRECTION_0") : c
        },
        _getStreetName: function(a) {
            var b = a.match(/\{highway:(.*)\}/);
            return b = b ? OSRM.loc("HIGHWAY_" + b[1].toUpperCase(), "HIGHWAY_DEFAULT") : a
        }
    }, OSRM.CONSTANTS.PRECISION = 6, OSRM.RoutingGeometry = {
        show: function(a) {
            var b = OSRM.RoutingGeometry._decode(a.route_geometry, OSRM.C.PRECISION);
            OSRM.G.route.showRoute(b, OSRM.Route.ROUTE)
        },
        showNA: function() {
            for (var a = [], b = 0, c = OSRM.G.markers.route.length; c > b; b++) a.push(OSRM.G.markers.route[b].getPosition());
            OSRM.G.route.showRoute(a, OSRM.Route.NOROUTE)
        },
        _decode: function(a, b) {
            b = Math.pow(10, -b);
            for (var c = a.length, d = 0, e = 0, f = 0, g = []; c > d;) {
                var h, i = 0,
                    j = 0;
                do h = a.charCodeAt(d++) - 63, j |= (31 & h) << i, i += 5; while (h >= 32);
                var k = 1 & j ? ~(j >> 1) : j >> 1;
                e += k, i = 0, j = 0;
                do h = a.charCodeAt(d++) - 63, j |= (31 & h) << i, i += 5; while (h >= 32);
                var l = 1 & j ? ~(j >> 1) : j >> 1;
                f += l, g.push([e * b, f * b])
            }
            return g
        }
    }, OSRM.GUI.extend({
        init: function() {
            OSRM.GUI.setDistanceFormat(OSRM.DEFAULTS.DISTANCE_FORMAT), document.getElementById("gui-input-source").onchange = function() {
                OSRM.GUI.inputChanged(OSRM.C.SOURCE_LABEL)
            }, document.getElementById("gui-delete-source").onclick = function() {
                OSRM.GUI.deleteMarker(OSRM.C.SOURCE_LABEL)
            }, document.getElementById("gui-search-source").onclick = function() {
                OSRM.GUI.showMarker(OSRM.C.SOURCE_LABEL)
            }, document.getElementById("gui-input-target").onchange = function() {
                OSRM.GUI.inputChanged(OSRM.C.TARGET_LABEL)
            }, document.getElementById("gui-delete-target").onclick = function() {
                OSRM.GUI.deleteMarker(OSRM.C.TARGET_LABEL)
            }, document.getElementById("gui-search-target").onclick = function() {
                OSRM.GUI.showMarker(OSRM.C.TARGET_LABEL)
            }, document.getElementById("gui-reset").onclick = OSRM.GUI.resetRouting, document.getElementById("gui-reverse").onclick = OSRM.GUI.reverseRouting, document.getElementById("open-editor").onclick = OSRM.GUI.openEditor, document.getElementById("open-josm").onclick = OSRM.GUI.openJOSM, document.getElementById("open-osmbugs").onclick = OSRM.GUI.openOSMBugs, document.getElementById("option-highlight-nonames").onclick = OSRM.GUI.hightlightNonames, document.getElementById("option-show-previous-routes").onclick = OSRM.GUI.showPreviousRoutes
        },
        activateRouteFeatures: function() {
            OSRM.Printing.activate(), OSRM.G.map.locationsControl.activateRoute(), OSRM.G.active_shortlink = null
        },
        deactivateRouteFeatures: function() {
            OSRM.Printing.deactivate(), OSRM.G.map.locationsControl.deactivateRoute(), OSRM.G.active_shortlink = null
        },
        resetRouting: function() {
            document.getElementById("gui-input-source").value = "", document.getElementById("gui-input-target").value = "", OSRM.G.route.reset(), OSRM.G.markers.reset(), OSRM.GUI.clearResults(), OSRM.JSONP.reset()
        },
        reverseRouting: function() {
            var a = document.getElementById("gui-input-source").value;
            document.getElementById("gui-input-source").value = document.getElementById("gui-input-target").value, document.getElementById("gui-input-target").value = a, OSRM.G.route.isShown() ? (OSRM.G.markers.route.reverse(), OSRM.Routing._Reversed(), OSRM.G.markers.route.reverse(), OSRM.G.markers.highlight.hide(), OSRM.RoutingDescription.showSimple(OSRM.G.response)) : OSRM.G.markers.reverseMarkers(), OSRM.G.markers.reverseDescriptions()
        },
        showMarker: function(a) {
            OSRM.JSONP.fences.geocoder_source || OSRM.JSONP.fences.geocoder_target || (a == OSRM.C.SOURCE_LABEL && OSRM.G.markers.hasSource() ? OSRM.G.markers.route[0].centerView() : a == OSRM.C.TARGET_LABEL && OSRM.G.markers.hasTarget() && OSRM.G.markers.route[OSRM.G.markers.route.length - 1].centerView())
        },
        inputChanged: function(a) {
            a == OSRM.C.SOURCE_LABEL ? OSRM.Geocoder.call(OSRM.C.SOURCE_LABEL, document.getElementById("gui-input-source").value) : a == OSRM.C.TARGET_LABEL && OSRM.Geocoder.call(OSRM.C.TARGET_LABEL, document.getElementById("gui-input-target").value)
        },
        openJOSM: function() {
            var a = OSRM.G.map.getZoom();
            if (a < OSRM.DEFAULTS.JOSM_MIN_ZOOM_LEVEL) window.alert(OSRM.loc("OPEN_JOSM_FAILED"));
            else {
                var b = OSRM.G.map.getCenterUI(),
                    c = OSRM.G.map.getBoundsUI(),
                    d = Math.min(.02, Math.abs(c.getSouthWest().lng - b.lng)),
                    e = Math.min(.01, Math.abs(c.getSouthWest().lat - b.lat)),
                    f = ["left=" + (b.lng - d).toFixed(6), "bottom=" + (b.lat - e).toFixed(6), "right=" + (b.lng + d).toFixed(6), "top=" + (b.lat + e).toFixed(6)],
                    g = "http://127.0.0.1:8111/load_and_zoom?" + f.join("&"),
                    h = document.getElementById("josm-frame");
                h || (h = L.DomUtil.create("iframe", null, document.body), h.style.display = "none", h.id = "josm-frame"), h.src = g
            }
        },
        openEditor: function() {
            var a = OSRM.G.map.getZoom();
            if (a < OSRM.DEFAULTS.EDITOR_MIN_ZOOM_LEVEL) window.alert(OSRM.loc("OPEN_EDITOR_FAILED"));
            else {
                var b = OSRM.G.map.getCenterUI(),
                    c = OSRM.C.PRECISION;
                window.open("http://www.openstreetmap.org/edit?lat=" + b.lat.toFixed(c) + "&lon=" + b.lng.toFixed(c) + "&zoom=" + a)
            }
        },
        openOSMBugs: function() {
            var a = OSRM.G.map.getZoom();
            if (a < OSRM.DEFAULTS.NOTES_MIN_ZOOM_LEVEL) window.alert(OSRM.loc("OPEN_OSMBUGS_FAILED"));
            else {
                var b = OSRM.G.map.getCenterUI(),
                    c = OSRM.C.PRECISION;
                window.open("http://www.openstreetmap.org/?lat=" + b.lat.toFixed(c) + "&lon=" + b.lng.toFixed(c) + "&zoom=" + a + "&notes=yes")
            }
        },
        deleteMarker: function(a) {
            var b = null;
            "source" == a && OSRM.G.markers.hasSource() ? b = 0 : "target" == a && OSRM.G.markers.hasTarget() && (b = OSRM.G.markers.route.length - 1), null != b && (OSRM.G.markers.removeMarker(b), OSRM.Routing.getRoute(), OSRM.G.markers.highlight.hide())
        },
        showPreviousRoutes: function() {
            0 == document.getElementById("option-show-previous-routes").checked ? OSRM.G.route.deactivateHistoryRoutes() : OSRM.G.route.activateHistoryRoutes()
        },
        zoomOnRoute: function() {
            if (0 != OSRM.G.route.isShown()) {
                var a = new L.LatLngBounds(OSRM.G.route._current_route.getPositions());
                OSRM.G.map.fitBoundsUI(a)
            }
        },
        zoomOnUser: function() {
            navigator.geolocation && navigator.geolocation.getCurrentPosition(OSRM.Map.geolocationResponse)
        },
        hightlightNonames: function() {
            OSRM.Routing.getRoute_Redraw({
                keepAlternative: !0
            })
        }
    }), OSRM.GUI.extend({
        init: function() {
            var a = OSRM.GUI.getRoutingEngines();
            OSRM.GUI.selectorInit("gui-engine-toggle", a, OSRM.DEFAULTS.ROUTING_ENGINE, OSRM.GUI._onRoutingEngineChanged)
        },
        setRoutingEngine: function(a) {
            a != OSRM.G.active_routing_engine && (OSRM.GUI.selectorChange("gui-engine-toggle", a), OSRM.G.active_routing_engine = a, OSRM.G.active_routing_metric = OSRM.DEFAULTS.ROUTING_ENGINES[OSRM.G.active_routing_engine].metric, OSRM.G.active_routing_server_url = OSRM.DEFAULTS.ROUTING_ENGINES[OSRM.G.active_routing_engine].url, OSRM.G.active_routing_timestamp_url = OSRM.DEFAULTS.ROUTING_ENGINES[OSRM.G.active_routing_engine].timestamp, OSRM.GUI.queryDataTimestamp())
        },
        _onRoutingEngineChanged: function(a) {
            a != OSRM.G.active_routing_engine && (OSRM.GUI.setRoutingEngine(a), OSRM.G.markers.route.length > 1 && OSRM.Routing.getRoute())
        },
        setRoutingEnginesLanguage: function() {
            var a = OSRM.GUI.getRoutingEngines();
            OSRM.GUI.selectorRenameOptions("gui-engine-toggle", a)
        },
        getRoutingEngines: function() {
            for (var a = OSRM.DEFAULTS.ROUTING_ENGINES, b = [], c = 0, d = a.length; d > c; c++) b.push({
                display: OSRM.loc(a[c].label),
                value: c
            });
            return b
        }
    }), OSRM.RoutingNoNames = {
        show: function(a) {
            if (0 == document.getElementById("option-highlight-nonames").checked) return void OSRM.G.route.hideUnnamedRoute();
            for (var b = [], c = 0; c < a.route_instructions.length; c++) b[a.route_instructions[c][3]] = "" == a.route_instructions[c][1] ? !1 : !0;
            for (var d = OSRM.RoutingGeometry._decode(a.route_geometry, OSRM.C.PRECISION), e = !0, f = [], g = [], c = 0; c < d.length; c++) f.push(d[c]), (b[c] != e && void 0 != b[c] || c == d.length - 1) && (0 == e && g.push(f), f = [], f.push(d[c]), e = b[c]);
            OSRM.G.route.showUnnamedRoute(g)
        },
        showNA: function() {
            OSRM.G.route.hideUnnamedRoute()
        }
    }, OSRM.Via = {
        _findNearestRouteSegment: function(a) {
            for (var b = Number.MAX_VALUE, c = void 0, d = OSRM.G.map.latLngToLayerPoint(a), e = OSRM.G.route.getPoints(), f = 1; f < e.length; f++) {
                var g = L.LineUtil._sqClosestPointOnSegment(d, e[f - 1], e[f], !0);
                b > g && (b = g, c = f)
            }
            return c
        },
        findViaIndex: function(a) {
            for (var b = OSRM.Via._findNearestRouteSegment(a), c = OSRM.G.response.via_points, d = c.length - 2, e = Array(), f = 1; f < c.length - 1; f++)
                if (e[f - 1] = OSRM.Via._findNearestRouteSegment(new L.LatLng(c[f][0], c[f][1])), e[f - 1] > b) {
                    d = f - 1;
                    break
                }
            return d
        },
        dragTimer: new Date,
        drawDragMarker: function(a) {
            if (0 != OSRM.G.route.isShown() && 1 != OSRM.G.dragging && !(new Date - OSRM.Via.dragTimer < 25)) {
                OSRM.Via.dragTimer = new Date;
                for (var b = OSRM.G.route._current_route.route.closestLayerPoint(a.layerPoint), c = b ? b.distance : 1e3, d = a.latlng, e = 0, f = OSRM.G.markers.route.length; f > e; e++)
                    if ("drag" != OSRM.G.markers.route[e].label) {
                        var g = OSRM.G.markers.route[e].getPosition(),
                            h = OSRM.G.map.project(d).distanceTo(OSRM.G.map.project(g));
                        20 > h && (c = 1e3)
                    }
                for (var i = OSRM.G.map.layerPointToContainerPoint(a.layerPoint), j = document.elementFromPoint(i.x, i.y), e = 0, f = OSRM.G.markers.route.length; f > e; e++) "drag" != OSRM.G.markers.route[e].label && j == OSRM.G.markers.route[e].marker._icon && (c = 1e3);
                OSRM.G.markers.highlight.isShown() && (OSRM.G.map.project(d).distanceTo(OSRM.G.map.project(OSRM.G.markers.highlight.getPosition())) < 20 ? c = 1e3 : j == OSRM.G.markers.highlight.marker._icon && (c = 1e3)), 20 > c ? (OSRM.G.markers.dragger.setPosition(OSRM.G.map.layerPointToLatLng(b)), OSRM.G.markers.dragger.show()) : OSRM.G.markers.dragger.hide()
            }
        }
    }, OSRM.CONSTANTS.SOURCE_LABEL = "source", OSRM.CONSTANTS.TARGET_LABEL = "target", OSRM.CONSTANTS.VIA_LABEL = "via", OSRM.CONSTANTS.DO_FALLBACK_TO_LAT_LNG = !0, OSRM.Geocoder = {
        call: function(a, b) {
            if ("" != b) {
                if (b.match(/^\s*[-+]?[0-9]*\.?[0-9]+\s*[,;]\s*[-+]?[0-9]*\.?[0-9]+\s*$/)) {
                    var c = b.split(/[,;]/);
                    return OSRM.Geocoder._onclickResult(a, c[0], c[1]), void OSRM.Geocoder.updateAddress(a)
                }
                var d = OSRM.DEFAULTS.HOST_GEOCODER_URL + "?format=json&json_callback=%jsonp" + OSRM.DEFAULTS.GEOCODER_BOUNDS + "&accept-language=" + OSRM.Localization.current_language + "&limit=30&q=" + b,
                    e = OSRM.G.map.getBounds();
                d += "&viewbox=" + e._southWest.lng + "," + e._northEast.lat + "," + e._northEast.lng + "," + e._southWest.lat, OSRM.JSONP.call(d, OSRM.Geocoder._showResults, OSRM.Geocoder._showResults_Timeout, OSRM.DEFAULTS.JSONP_TIMEOUT, "geocoder_" + a, {
                    marker_id: a,
                    query: b
                })
            }
        },
        _onclickResult: function(a, b, c, d) {
            var e;
            if (a == OSRM.C.SOURCE_LABEL) e = OSRM.G.markers.setSource(new L.LatLng(b, c));
            else {
                if (a != OSRM.C.TARGET_LABEL) return;
                e = OSRM.G.markers.setTarget(new L.LatLng(b, c))
            }
            OSRM.G.markers.route[e].show(), OSRM.G.markers.route[e].centerView(d), OSRM.G.markers.route.length > 1 && OSRM.Routing.getRoute()
        },
        _showResults: function(a, b) {
            if (!a) return void OSRM.Geocoder._showResults_Empty(b);
            if (0 == a.length) return void OSRM.Geocoder._showResults_Empty(b);
            var c = a,
                d = null;
            if (null != c[0].boundingbox && (d = OSRM.G.map.getBoundsZoom([
                    [c[0].boundingbox[0], c[0].boundingbox[2]],
                    [c[0].boundingbox[1], c[0].boundingbox[3]]
                ], !0)), OSRM.Geocoder._onclickResult(b.marker_id, c[0].lat, c[0].lon, d), !(OSRM.G.markers.route.length > 1)) {
                var e = "";
                e += '<table class="results medium-font">';
                for (var f = 0; f < c.length; f++) {
                    var g = c[f],
                        h = "results-body-odd";
                    if (f % 2 == 0 && (h = "results-body-even"), e += '<tr class="' + h + '">', g.icon || (g.icon = "http://nominatim.openstreetmap.org/images/mapicons/poi_point_of_interest.glow.12.png"), e += '<td class="results-body-counter"><img src="' + g.icon + '" alt=""/></td>', e += '<td class="results-body-items">', g.display_name) {
                        var d = "null";
                        null != g.boundingbox && (d = OSRM.G.map.getBoundsZoom([
                            [g.boundingbox[0], g.boundingbox[2]],
                            [g.boundingbox[1], g.boundingbox[3]]
                        ], !0)), e += '<div class="results-body-item" onclick="OSRM.Geocoder._onclickResult(\'' + b.marker_id + "', " + g.lat + ", " + g.lon + ", " + d + ');">' + g.display_name, e += "</div>"
                    }
                    e += "</td></tr>"
                }
                e += "</table>", document.getElementById("information-box-header").innerHTML = "<div class='header-title'>" + OSRM.loc("SEARCH_RESULTS") + "</div><div class='header-content'>(" + OSRM.loc("FOUND_X_RESULTS").replace(/%i/, c.length) + ")</div>", "<div class='header-content'>(found " + c.length + " results)</div>", document.getElementById("information-box").className = "information-box-with-normal-header", document.getElementById("information-box").innerHTML = e
            }
        },
        _showResults_Empty: function(a) {
            document.getElementById("information-box-header").innerHTML = "<div class='header-title'>" + OSRM.loc("SEARCH_RESULTS") + "</div><div class='header-content'>(" + OSRM.loc("FOUND_X_RESULTS").replace(/%i/, 0) + ")</div>", document.getElementById("information-box").className = "information-box-with-normal-header", document.getElementById("information-box").innerHTML = a.marker_id == OSRM.C.SOURCE_LABEL ? "<div class='no-results big-font'>" + OSRM.loc("NO_RESULTS_FOUND_SOURCE") + ": " + a.query + "</div>" : a.marker_id == OSRM.C.TARGET_LABEL ? "<div class='no-results big-font'>" + OSRM.loc("NO_RESULTS_FOUND_TARGET") + ": " + a.query + "</div>" : "<div class='no-results big-font'>" + OSRM.loc("NO_RESULTS_FOUND") + ": " + a.query + "</div>"
        },
        _showResults_Timeout: function() {
            document.getElementById("information-box-header").innerHTML = "<div class='header-title'>" + OSRM.loc("SEARCH_RESULTS") + "</div><div class='header-content'>(" + OSRM.loc("FOUND_X_RESULTS").replace(/%i/, 0) + ")</div>", document.getElementById("information-box").className = "information-box-with-normal-header", document.getElementById("information-box").innerHTML = "<div class='no-results big-font'>" + OSRM.loc("TIMED_OUT") + "</div>"
        },
        updateLocation: function(a) {
            var b = OSRM.C.PRECISION;
            a == OSRM.C.SOURCE_LABEL && OSRM.G.markers.hasSource() ? document.getElementById("gui-input-source").value = OSRM.G.markers.route[0].getLat().toFixed(b) + ", " + OSRM.G.markers.route[0].getLng().toFixed(b) : a == OSRM.C.TARGET_LABEL && OSRM.G.markers.hasTarget() && (document.getElementById("gui-input-target").value = OSRM.G.markers.route[OSRM.G.markers.route.length - 1].getLat().toFixed(b) + ", " + OSRM.G.markers.route[OSRM.G.markers.route.length - 1].getLng().toFixed(b))
        },
        updateAddress: function(a, b) {
            var c = null,
                d = null,
                e = null;
            if (a == OSRM.C.SOURCE_LABEL && OSRM.G.markers.hasSource()) c = OSRM.G.markers.route[0].getLat(), d = OSRM.G.markers.route[0].getLng(), e = OSRM.G.markers.route[0].description;
            else {
                if (a != OSRM.C.TARGET_LABEL || !OSRM.G.markers.hasTarget()) return;
                c = OSRM.G.markers.route[OSRM.G.markers.route.length - 1].getLat(), d = OSRM.G.markers.route[OSRM.G.markers.route.length - 1].getLng(), e = OSRM.G.markers.route[OSRM.G.markers.route.length - 1].description
            }
            if (null != e) return void OSRM.Geocoder._showReverseResults({
                address: {
                    road: e
                }
            }, {
                marker_id: a
            });
            var f = OSRM.C.PRECISION,
                g = OSRM.DEFAULTS.HOST_REVERSE_GEOCODER_URL + "?format=json&json_callback=%jsonp&accept-language=" + OSRM.Localization.current_language + "&lat=" + c.toFixed(f) + "&lon=" + d.toFixed(f);
            OSRM.JSONP.call(g, OSRM.Geocoder._showReverseResults, OSRM.Geocoder._showReverseResults_Timeout, OSRM.DEFAULTS.JSONP_TIMEOUT, "reverse_geocoder_" + a, {
                marker_id: a,
                do_fallback: b
            })
        },
        _showReverseResults: function(a, b) {
            if (!a) return void OSRM.Geocoder._showReverseResults_Timeout(a, b);
            if (void 0 == a.address) return void OSRM.Geocoder._showReverseResults_Timeout(a, b);
            var c = 0,
                d = "";
            return a.address.road && (d += a.address.road, c++), a.address.city ? (c > 0 && (d += ", "), d += a.address.city, c++) : a.address.village && (c > 0 && (d += ", "), d += a.address.village, c++), 2 > c && a.address.country && (c > 0 && (d += ", "), d += a.address.country, c++), 0 == c ? void OSRM.Geocoder._showReverseResults_Timeout(a, b) : void(b.marker_id == OSRM.C.SOURCE_LABEL && OSRM.G.markers.hasSource() ? document.getElementById("gui-input-source").value = d : b.marker_id == OSRM.C.TARGET_LABEL && OSRM.G.markers.hasTarget() && (document.getElementById("gui-input-target").value = d))
        },
        _showReverseResults_Timeout: function(a, b) {
            b.do_fallback && OSRM.Geocoder.updateLocation(b.marker_id)
        },
        _showInitResults: function(a, b) {
            var c = OSRM.G.initial_positions;
            if (c.done++, a && 0 != a.length ? (c.positions[b.id] = new L.LatLng(a[0].lat, a[0].lon), null == c.zoom && null != a[0].boundingbox && (c.zoom = OSRM.G.map.getBoundsZoom([
                    [a[0].boundingbox[0], a[0].boundingbox[2]],
                    [a[0].boundingbox[1], a[0].boundingbox[3]]
                ], !0))) : c.fail = !0, c.done == c.positions.length) {
                if (1 == c.fail) return void OSRM.GUI.exclusiveNotify(OSRM.loc("NOTIFICATION_GEOCODERFAIL_HEADER"), OSRM.loc("NOTIFICATION_GEOCODERFAIL_BODY"), !0);
                OSRM.GUI.exclusiveDenotify(), OSRM.Geocoder[b.callback]()
            }
        },
        _showInitResults_Destinations: function() {
            var a = OSRM.G.initial_positions,
                b = a.positions,
                c = OSRM.G.markers.setTarget(b[b.length - 1]);
            OSRM.G.markers.route[c].description = a.name, OSRM.Geocoder.updateAddress(OSRM.C.TARGET_LABEL, OSRM.C.DO_FALLBACK_TO_LAT_LNG), OSRM.G.markers.route[c].show(), OSRM.G.markers.route[c].centerView(a.zoom);
            for (var d = 0; d < b.length - 1; d++) OSRM.G.markers.addInitialVia(b[d]);
            OSRM.G.initial_position_override = !0
        },
        _showInitResults_Locations: function() {
            var a = OSRM.G.initial_positions,
                b = a.positions;
            b.length > 0 && (OSRM.G.markers.setSource(b[0]), OSRM.Geocoder.updateAddress(OSRM.C.SOURCE_LABEL, OSRM.C.DO_FALLBACK_TO_LAT_LNG)), b.length > 1 && (OSRM.G.markers.setTarget(b[b.length - 1]), OSRM.Geocoder.updateAddress(OSRM.C.TARGET_LABEL, OSRM.C.DO_FALLBACK_TO_LAT_LNG));
            for (var c = 1; c < b.length - 1; c++) OSRM.G.markers.setVia(c - 1, b[c]);
            for (var c = 0; c < OSRM.G.markers.route.length; c++) OSRM.G.markers.route[c].show();
            if (null == a.zoom || null == a.center) {
                var d = new L.LatLngBounds(b);
                OSRM.G.map.fitBoundsUI(d)
            } else OSRM.G.map.setView(a.center, a.zoom);
            OSRM.G.active_alternative = a.active_alternative || 0, OSRM.Routing.getRoute({
                keepAlternative: !0
            }), OSRM.G.initial_position_override = !0
        }
    }, OSRM.CSS = {
        getStylesheet: function(a, b) {
            b = b || document;
            for (var c = b.styleSheets, d = 0, e = c.length; e > d; d++)
                if (c[d].href.indexOf(a) >= 0) return c[d];
            return null
        },
        insert: function(a, b, c) {
            a.addRule ? a.addRule(b, c) : a.insertRule && a.insertRule(b + " { " + c + " }", a.cssRules.length)
        }
    }, OSRM.JSONP = {
        fences: {},
        callbacks: {},
        timeouts: {},
        timers: {},
        late: function() {},
        empty: function() {},
        call: function(a, b, c, d, e, f) {
            if (1 == OSRM.JSONP.fences[e]) return !1;
            OSRM.JSONP.fences[e] = !0, OSRM.JSONP.timeouts[e] = function(a) {
                try {
                    c(a, f)
                } finally {
                    OSRM.JSONP.callbacks[e] = OSRM.JSONP.late, OSRM.JSONP.timeouts[e] = OSRM.JSONP.empty, OSRM.JSONP.fences[e] = void 0
                }
            }, OSRM.JSONP.callbacks[e] = function(a) {
                clearTimeout(OSRM.JSONP.timers[e]), OSRM.JSONP.timers[e] = void 0;
                try {
                    b(a, f)
                } finally {
                    OSRM.JSONP.callbacks[e] = OSRM.JSONP.empty, OSRM.JSONP.timeouts[e] = OSRM.JSONP.late, OSRM.JSONP.fences[e] = void 0
                }
            };
            var g = document.getElementById("jsonp_" + e);
            g && g.parentNode.removeChild(g);
            var h = document.createElement("script");
            return h.type = "text/javascript", h.id = "jsonp_" + e, h.src = a.replace(/%jsonp/, "OSRM.JSONP.callbacks." + e), document.head.appendChild(h), OSRM.JSONP.timers[e] = setTimeout(OSRM.JSONP.timeouts[e], d), !0
        },
        clear: function(a) {
            clearTimeout(OSRM.JSONP.timers[a]), OSRM.JSONP.callbacks[a] = OSRM.JSONP.empty, OSRM.JSONP.timeouts[a] = OSRM.JSONP.empty, OSRM.JSONP.fences[a] = void 0;
            var b = document.getElementById("jsonp_" + a);
            b && b.parentNode.removeChild(b)
        },
        reset: function() {
            OSRM.JSONP.fences = {}, OSRM.JSONP.callbacks = {}, OSRM.JSONP.timeouts = {}, OSRM.JSONP.timers = {}
        }
    }, OSRM.Localization = {
        DIRECTORY: "localization/",
        current_language: OSRM.DEFAULTS.LANGUAGE,
        fallback_language: "en",
        load_on_demand_language: null,
        init: function() {
            var a = OSRM.DEFAULTS.LANGUAGE_SUPPORTED;
            if (1 == OSRM.DEFAULTS.LANGUAGE_USE_BROWSER_SETTING)
                for (var b = (navigator.language || navigator.userLanguage || "").substring(0, 2), c = 0; c < a.length; ++c) a[c].encoding == b && (OSRM.Localization.current_language = b);
            for (var d = [], e = [], f = 0, c = 0, g = a.length; g > c; c++) d.push({
                display: a[c].encoding,
                value: a[c].encoding
            }), e.push({
                display: a[c].name,
                value: a[c].encoding
            }), a[c].encoding == OSRM.Localization.current_language && (f = c);
            OSRM.GUI.selectorInit("gui-language-toggle", d, f, OSRM.Localization.setLanguageWrapper), OSRM.GUI.selectorInit("gui-language-2-toggle", e, f, OSRM.Localization.setLanguageWrapper), OSRM.Localization.setLanguage(OSRM.Localization.fallback_language), OSRM.Localization.setLanguage(OSRM.Localization.current_language)
        },
        setLanguageWrapper: function(a) {
            OSRM.GUI.deactivateTooltip("LOCALIZATION"), OSRM.Localization.setLanguage(a)
        },
        setLanguage: function(a, b) {
            if (b) {
                var c = OSRM.Localization[OSRM.Localization.fallback_language];
                if (null == c || null != c.loading) return;
                var d = OSRM.Localization[OSRM.Localization.load_on_demand_language];
                if (null != d && null == d.loading && (a = OSRM.Localization.load_on_demand_language), a != OSRM.Localization.load_on_demand_language) return
            }
            if (OSRM.GUI.selectorChange("gui-language-toggle", a), OSRM.GUI.selectorChange("gui-language-2-toggle", a), OSRM.Localization[a]) {
                if (OSRM.Localization[a].loading) return;
                OSRM.Localization.current_language = a, OSRM.Localization.load_on_demand_language = null, OSRM.GUI.setLabels(), OSRM.GUI.updateNotifications(), OSRM.Utils.updateAbbreviationCache();
                for (var e = 0, f = OSRM.G.localizable_maps.length; f > e; e++) OSRM.G.localizable_maps[e].options.culture = OSRM.loc("CULTURE");
                if (OSRM.G.map.layerControl.getActiveLayer().redraw && OSRM.G.map.layerControl.getActiveLayer().redraw(), OSRM.G.map.layerControl.setLayerLabels(), null == OSRM.G.markers) return;
                OSRM.G.markers.route.length > 1 ? OSRM.Routing.getRoute({
                    keepAlternative: !0
                }) : OSRM.G.markers.route.length > 0 && "" != document.getElementById("information-box").innerHTML ? (OSRM.Geocoder.call(OSRM.C.SOURCE_LABEL, document.getElementById("gui-input-source").value), OSRM.Geocoder.call(OSRM.C.TARGET_LABEL, document.getElementById("gui-input-target").value)) : (OSRM.Geocoder.updateAddress(OSRM.C.SOURCE_LABEL, !1), OSRM.Geocoder.updateAddress(OSRM.C.TARGET_LABEL, !1), OSRM.GUI.clearResults())
            } else if (1 == OSRM.DEFAULTS.LANUGAGE_ONDEMAND_RELOADING)
                for (var g = OSRM.DEFAULTS.LANGUAGE_SUPPORTED, e = 0, f = g.length; f > e; e++)
                    if (g[e].encoding == a) {
                        OSRM.Localization.load_on_demand_language = a, OSRM.Localization[a] = {
                            loading: !0
                        };
                        var h = document.createElement("script");
                        h.type = "text/javascript", h.src = OSRM.Localization.DIRECTORY + "OSRM.Locale." + a + ".js", document.head.appendChild(h);
                        break
                    }
        },
        translate: function(a, b) {
            return OSRM.Localization[OSRM.Localization.current_language] && OSRM.Localization[OSRM.Localization.current_language][a] ? OSRM.Localization[OSRM.Localization.current_language][a] : OSRM.Localization[OSRM.Localization.fallback_language] && OSRM.Localization[OSRM.Localization.fallback_language][a] ? OSRM.Localization[OSRM.Localization.fallback_language][a] : b ? OSRM.loc(b) : a
        }
    }, OSRM.loc = OSRM.Localization.translate, OSRM.Printing = {
        DIRECTORY: "printing/",
        BASE_DIRECTORY: "../",
        init: function() {
            var a = document.createElement("div");
            a.id = "gui-printer-inactive", a.className = "iconic-button top-right-button";
            var b = document.createElement("div");
            b.className = "quad top-right-button";
            var c = document.getElementById("input-mask-header");
            c.appendChild(b, c.lastChild), c.appendChild(a, c.lastChild), document.getElementById("gui-printer-inactive").onclick = OSRM.Printing.openPrintWindow, OSRM.Browser.onUnloadHandler(OSRM.Printing.uninit)
        },
        uninit: function() {
            OSRM.G.printwindow && OSRM.G.printwindow.close()
        },
        activate: function() {
            document.getElementById("gui-printer-inactive") && (document.getElementById("gui-printer-inactive").id = "gui-printer")
        },
        deactivate: function() {
            document.getElementById("gui-printer") && (document.getElementById("gui-printer").id = "gui-printer-inactive")
        },
        show: function(a) {
            for (var b = "(", c = 0, d = a.route_name.length; d > c; c++) b += (c > 0 && "" != a.route_name[c] && "" != a.route_name[c - 1] ? " - " : "") + "<span style='white-space:nowrap;'>" + a.route_name[c] + "</span>";
            b += ")";
            var e;
            e = OSRM.Browser.IE6_8 ? '<thead class="description-header"><tr><td colspan="3"><table class="full"><tr class="row"><td class="left stretch"><table class="full"><tr class="row"><td class="center description-header-label nowrap">' + OSRM.loc("GUI_START") + ': </td><td class="left description-header-content stretch">' + document.getElementById("gui-input-source").value + '</td></tr><tr class="row"><td class="center description-header-label nowrap">↓</td><td class="left description-header-label stretch">' + b + '</td></tr><tr class="row"><td class="center description-header-label nowrap">' + OSRM.loc("GUI_END") + ': </td><td class="left description-header-content stretch">' + document.getElementById("gui-input-target").value + '</td></tr></table></td><td class="left"><table class="full"><tr class="row"><td class="left description-header-label nowrap">' + OSRM.loc("DISTANCE") + ': </td><td class="left description-header-content stretch">' + OSRM.Utils.toHumanDistance(a.route_summary.total_distance) + '</td></tr><tr class="row"><td class="left description-header-label nowrap">' + OSRM.loc("DURATION") + ': </td><td class="left description-header-content stretch">' + OSRM.Utils.toHumanTime(a.route_summary.total_time) + '</td></tr></table></td></tr></table><div class="quad"></div></td></tr></thead>' : '<thead class="description-header"><tr><td colspan="3"><div class="full"><div class="row"><div class="left stretch"><div class="full"><div class="row"><div class="center description-header-label nowrap">' + OSRM.loc("GUI_START") + ': </div><div class="left description-header-content stretch">' + document.getElementById("gui-input-source").value + '</div></div><div class="row"><div class="center description-header-label nowrap">↓</div><div class="left description-header-label stretch">' + b + '</div></div><div class="row"><div class="center description-header-label nowrap">' + OSRM.loc("GUI_END") + ': </div><div class="left description-header-content stretch">' + document.getElementById("gui-input-target").value + '</div></div></div></div><div class="left"><div class="full"><div class="row"><div class="left description-header-label nowrap">' + OSRM.loc("DISTANCE") + ': </div><div class="left description-header-content stretch">' + OSRM.Utils.toHumanDistance(a.route_summary.total_distance) + '</div></div><div class="row"><div class="left description-header-label nowrap">' + OSRM.loc("DURATION") + ': </div><div class="left description-header-content stretch">' + OSRM.Utils.toHumanTime(a.route_summary.total_time) + '</div></div></div></div></div></div><div class="quad"></div></td></tr></thead>';
            for (var f = '<tbody class="description-body">', g = 0; g < a.route_instructions.length; g++) {
                var h = "description-body-odd";
                g % 2 == 0 && (h = "description-body-even"), f += '<tr class="' + h + '">', f += '<td class="description-body-directions">', f += '<img class="description-body-direction" src="' + OSRM.Printing.BASE_DIRECTORY + OSRM.RoutingDescription._getDrivingInstructionIcon(a.route_instructions[g][0]) + '" alt="" />', f += "</td>", f += '<td class="description-body-items">', f += "" != a.route_instructions[g][1] ? OSRM.loc(OSRM.RoutingDescription._getDrivingInstruction(a.route_instructions[g][0])).replace(/\[(.*)\]/, "$1").replace(/%s/, a.route_instructions[g][1]).replace(/%d/, OSRM.loc(a.route_instructions[g][6])) : OSRM.loc(OSRM.RoutingDescription._getDrivingInstruction(a.route_instructions[g][0])).replace(/\[(.*)\]/, "").replace(/%d/, OSRM.loc(a.route_instructions[g][6])), f += "</td>", f += '<td class="description-body-distance">', f += g == a.route_instructions.length - 1 ? "&nbsp;" : "<b>" + OSRM.Utils.toHumanDistance(a.route_instructions[g][2]) + "</b>", f += "</td>", f += "</tr>"
            }
            f += "</tbody>";
            var i = OSRM.G.printwindow;
            i.document.getElementById("description").innerHTML = '<table class="description medium-font">' + e + f + "</table>", i.document.getElementById("overview-map-description").innerHTML = '<table class="description medium-font">' + e + "</table>";
            var j = OSRM.G.map.getActiveLayerId(),
                k = OSRM.G.route.getPositions(),
                l = new L.LatLngBounds(k),
                m = i.OSRM.drawMap(OSRM.DEFAULTS.TILE_SERVERS[j], l);
            i.OSRM.prefetchIcons(OSRM.G.images), i.OSRM.drawMarkers(OSRM.G.markers.route), i.OSRM.drawRoute(k), OSRM.JSONP.call(OSRM.Routing._buildCall() + "&z=" + m + "&instructions=false", OSRM.Printing.drawRoute, OSRM.Printing.timeoutRoute, OSRM.DEFAULTS.JSONP_TIMEOUT, "print")
        },
        timeoutRoute: function() {},
        drawRoute: function(a) {
            a && (a.alternative_geometries.unshift(a.route_geometry), OSRM.G.active_alternative >= a.alternative_geometries.length || (positions = OSRM.RoutingGeometry._decode(a.alternative_geometries[OSRM.G.active_alternative], OSRM.C.PRECISION), OSRM.G.printwindow.OSRM.drawRoute(positions)))
        },
        openPrintWindow: function() {
            OSRM.G.route.isRoute() && OSRM.G.route.isShown() && (OSRM.G.printwindow && OSRM.G.printwindow.close(), OSRM.G.printwindow = window.open(OSRM.Printing.DIRECTORY + "printing.html", "", "width=540,height=500,left=100,top=100,dependent=yes,location=no,menubar=no,scrollbars=yes,status=no,toolbar=no,resizable=yes"))
        },
        printWindowLoaded: function() {
            for (var a = OSRM.G.printwindow, b = a.document, c = [{
                    id: "#gui-printer-inactive",
                    image_id: "printer_inactive"
                }, {
                    id: "#gui-printer",
                    image_id: "printer"
                }, {
                    id: "#gui-printer:hover",
                    image_id: "printer_hover"
                }, {
                    id: "#gui-printer:active",
                    image_id: "printer_active"
                }], d = OSRM.CSS.getStylesheet("printing.css", b), e = 0; e < c.length; e++) OSRM.CSS.insert(d, c[e].id, "background-image:url(" + OSRM.Printing.BASE_DIRECTORY + OSRM.G.images[c[e].image_id].getAttribute("src") + ");");
            return a.OSRM.G.active_distance_format = OSRM.G.active_distance_format, a.OSRM.G.Localization.culture = OSRM.loc("CULTURE"), b.getElementById("description-label").innerHTML = OSRM.loc("ROUTE_DESCRIPTION"), b.getElementById("overview-map-label").innerHTML = OSRM.loc("OVERVIEW_MAP"), OSRM.G.route.isRoute() && OSRM.G.route.isShown() ? (OSRM.Printing.show(OSRM.G.response), b.getElementById("gui-printer-inactive").id = "gui-printer", b.getElementById("gui-printer").onclick = a.printWindow, void a.focus()) : void(b.getElementById("description").innerHTML = OSRM.loc("NO_ROUTE_SELECTED"))
        }
    }, OSRM.Utils = {
        seconds: "s",
        minutes: "min",
        hours: "h",
        miles: "mi",
        feet: "ft",
        kilometers: "km",
        meters: "m",
        updateAbbreviationCache: function() {
            OSRM.Utils.seconds = OSRM.loc("GUI_S"), OSRM.Utils.minutes = OSRM.loc("GUI_MIN"), OSRM.Utils.hours = OSRM.loc("GUI_H"), OSRM.Utils.miles = OSRM.loc("GUI_MI"), OSRM.Utils.feet = OSRM.loc("GUI_FT"), OSRM.Utils.kilometers = OSRM.loc("GUI_KM"), OSRM.Utils.meters = OSRM.loc("GUI_M")
        },
        toHumanTime: function(a) {
            return a = parseInt(a), minutes = parseInt(a / 60), a %= 60, hours = parseInt(minutes / 60), minutes %= 60, 0 == hours && 0 == minutes ? a + "&nbsp;s" : 0 == hours ? minutes + "&nbsp;min" : hours > 0 ? hours + "&nbsp;h&nbsp;" + minutes + "&nbsp;min" : "N/A"
        },
        toHumanDistanceMeters: function(a) {
            var b = parseInt(a);
            return b /= 1e3, b >= 100 ? b.toFixed(0) + "&nbsp;" + OSRM.Utils.kilometers : b >= 10 ? b.toFixed(1) + "&nbsp;" + OSRM.Utils.kilometers : b >= .1 ? b.toFixed(2) + "&nbsp;" + OSRM.Utils.kilometers : b >= 0 ? (1e3 * b).toFixed(0) + "&nbsp;" + OSRM.Utils.meters : "N/A"
        },
        toHumanDistanceMiles: function(a) {
            var b = parseInt(a);
            return b /= 1609.344, b >= 100 ? b.toFixed(0) + "&nbsp;" + OSRM.Utils.miles : b >= 10 ? b.toFixed(1) + "&nbsp;" + OSRM.Utils.miles : b >= .1 ? b.toFixed(2) + "&nbsp;" + OSRM.Utils.miles : b >= 0 ? (5280 * b).toFixed(0) + "&nbsp;" + OSRM.Utils.feet : "N/A"
        },
        toHumanDistance: null,
        isLatitude: function(a) {
            return a >= -90 && 90 >= a ? !0 : !1
        },
        isLongitude: function(a) {
            return a >= -180 && 180 >= a ? !0 : !1
        },
        isNumber: function(a) {
            return !isNaN(parseFloat(a)) && isFinite(a)
        },
        getHash: function(a) {
            return a.split("").reduce(function(a, b) {
                return a = (a << 5) - a + b.charCodeAt(0), a & a
            }, 0)
        }
    };

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
// INT_MAX=2147483647;

// function compare_three_tree(a,b,c){
	// ok_count=0;
	// a_count=0;
	// b_count=0;
	// for(var i in glob.forward[c]) {
		// if(glob.backward[a][i]) {
			// ++a_count;
			// if(glob.backward[b][i]) {
				// ++ok_count;
				// ++b_count;
			// }
		// }
		// else if(glob.backward[b][i])
			// ++b_count;
	// }
	// return max(a_count, b_count) - ok_count;
// }
// function compare_tree(a,b){
	// ok_count=0;
	// str_ok_count=0;
	// min_dist=INT_MAX;
	// for(var i in glob.backward[b]) {
		// if(glob.backward[a][i]) ++ok_count;
		// if(glob.forward[a][i]){
			// ++str_ok_count;
			// if(min_dist>glob.forward[a][i].distance+glob.backward[b][i].distance)
				// min_dist=glob.forward[a][i].distance+glob.backward[b][i].distance;
		// }
	// }
	// return [ok_count, Math.abs(min_dist), Object.keys(glob.backward[a]).length, Object.keys(glob.backward[b]).length, str_ok_count];
// }
// function max(x1,x2) {
	// return x1>x2?x1:x2;
// }
// function balance(arr){
	// return  max(arr[2],arr[3])/arr[0]+P.balance*arr[1];
// }
// P={balance:1, nearest_count_radius:2, max_chains_from_point:INT_MAX,chain_balance:0.9};
// passing_map={};
// chains=[];
// points_heap=[];
// glob_pulls=[];
// priority=[];
// function add_chain_re(chain_pull, chain) {
	// var cur_node=chain[chain.length-1]
	// var from_cache=[];
	// if(glob_pulls[cur_node])
		// for(var i=0;i<glob_pulls[cur_node].length;++i) {
			// var ok=true;
			// for(var j=0; ok && j<chain.length-1; ++j)	
				// for(var k=0;k<glob_pulls[cur_node][i].length;++k)
					// if(chain[j]==glob_pulls[cur_node][i][k]) {
						// ok=false;
						// break;
					// }
			// if(ok)
				// for(var k in glob_pulls[cur_node][i].stopers) {
					// var stoper_in_chain=false;
					// for(var j=0; j<chain.length-1; ++j)
						// if(chain[j]==(k-0)) {
							// stoper_in_chain=true;
							// break;
						// }
					// if(!stoper_in_chain) {
						// ok=false;
						// break;
					// }
				// }
			// if(ok) 
				// from_cache.push(i);
		// }
	// if(from_cache.length) 
	// {
		// //for(var i=0;i<from_cache.length;++i){
		// if(from_cache.length>1) {
			// chain.cont=from_cache;
			// chain_pull.push(chain);
		// }
		// else{
			// var i=from_cache[0];
			// var new_chain=chain.concat(glob_pulls[cur_node][i]);
			// new_chain.stopers={};
			// for (var k in chain.stopers) { new_chain.stopers[k] = chain.stopers[k]; }
			// for (var k in glob_pulls[cur_node][i].stopers) { new_chain.stopers[k] = glob_pulls[cur_node][i].stopers[k]; }
			// new_chain.passing=chain.passing.concat(glob_pulls[cur_node][i].passing);
			// new_chain.passing_sum=chain.passing_sum+glob_pulls[cur_node][i].passing_sum;
			// chain_pull.push(new_chain);
		// }
		// return;
	// }
	// console.log(chain);
	// //var next_points=[];
	// var new_next_points=[];
	// //var cur_stopers=JSON.parse(JSON.stringify(chain.stopers))
	// for(var i in passing_map[cur_node])
		// if(chain.indexOf(i-0)==-1) new_next_points.push(i-0);
		// else chain.stopers[i]=true;
	// /*if(next_points.length<3) new_next_points=next_points;
	// else {
		// var mean=0;
		// for(var i=1; i<next_points.length; ++i)
			// mean+=compare_three_tree(next_points[0],next_points[i],cur_node)
		// mean/=next_points.length-1;
		// new_next_points.push(next_points[0]);
		// for(var i=1; i<next_points.length; ++i) {
			// var add=true;
			// for(var j=0; j<new_next_points.length; ++j)
				// if(compare_three_tree(new_next_points[j],next_points[i],cur_node)<mean) {
						// add=false;
						// if(passing_map[cur_node][new_next_points[j]]>passing_map[cur_node][next_points[i]]) {
							// new_next_points[j]=next_points[i];
							// break;
						// }
					// }
			// if(add)
				// new_next_points.push(next_points[i]);
		// }
	// }*/
		
	// if(new_next_points.length>0){
		// var cur_pull=[]
		// for(var i in new_next_points) {
			// var new_chain=chain.slice(0);
			// new_chain.push(new_next_points[i]);
			// new_chain.stopers=chain.stopers;
			// new_chain.passing=chain.passing.slice(0);
			// new_chain.passing.push(passing_map[cur_node][new_next_points[i]]);
			// new_chain.passing_sum=chain.passing_sum+passing_map[cur_node][new_next_points[i]];
			// add_chain_re(cur_pull, new_chain)
		// }
		// if(!glob_pulls[cur_node]) glob_pulls[cur_node]=[];
		// for(var i in cur_pull) {
			// chain_pull.push(cur_pull[i]);
			// var sub_chain=cur_pull[i].slice(chain.length);
			// sub_chain.stopers={};//JSON.parse(JSON.stringify(cur_pull[i].stopers));
			// for(var j=0;j<chain.length-1;++j)
				// if(cur_pull[i].stopers[chain[j]])
					// sub_chain.stopers[chain[j]]=true;
			// //for(var j=0;j<sub_chain.length;++j)
			// //		sub_chain.stopers[sub_chain[j]]=true;
			// sub_chain.passing=cur_pull[i].passing.slice(chain.passing.length);
			// sub_chain.passing_sum=cur_pull[i].passing_sum-chain.passing_sum;
			// glob_pulls[cur_node].push(sub_chain);
		// }
	// }
	// else
		// chain_pull.push(chain);
// }
// function clusterize(){

	// passing_map={};
	// //passing_map[glob.n]={};
	// //passing_map[glob.n][glob.n]=INT_MAX;
	// chains=[];
// glob_pulls=[];
	// for(var i=0;i<glob.n;++i) {
		// passing_map[i]={};
		// var temp=[]
		// for(var j=0;j<glob.n;++j)
			// if(i != j) 
				// temp.push([j, balance(compare_tree(i,j))]);
		// temp.sort(function(a,b){return a[1]-b[1]})
		// temp=temp.slice(0, P.nearest_count_radius);
		// for(var j=0;j<temp.length;++j)
			// passing_map[i][temp[j][0]]=temp[j][1];
	// }
	// for(var a=0; a<glob.n; ++a) {
		// var chain_pull=[];
		// var base_chain=[a];
		// base_chain.passing=[];
		// base_chain.passing_sum=0;
		// base_chain.stopers={};
		// add_chain_re(chain_pull, base_chain);
		// //chain_pull.sort(function(a,b) { return P.chain_balance*(a.length-b.length)+(1-P.chain_balance)*(a.passing_sum-b.passing_sum) });
		// chains=chains.concat(chain_pull.slice(0, P.max_chains_from_point));
	// }
	// showChains()
// }
function showChains(a){
	window.glob=a
	window.chains=window.glob.chains;
	var box=document.getElementById('clusters-box');
	var text='Выделить: <a href="#" onclick="return showAllMarkers()">все</a> <a href="#" onclick="return showMarkersByPriority()">значимые</a> <a href="#" onclick="return showMarkersByPriorityDesc()">редкие</a><ul>';
	for(var i=0;i<chains.length;i++) {
		text+='<li><a href="#" onclick="return showChain('+i+')">'+i+'</a>=';
		text+=chains[i].join('+');
		//text+=' ('+chains[i].passing_sum/chains[i].length+') '
		//if(chains[i].cont) text+=chains[i].cont.join();
		text+='</li>';
	}
	text+='</ul>';
	box.innerHTML=text;
}
function showChain(cur){
	for(var i=0; i<glob.n; ++i) {
		OSRM.G.markers.route[i].marker.setIcon(OSRM.G.icons['marker-highlight'])
		OSRM.G.markers.route[i].marker.setOpacity(0.8);
	}
	var b = OSRM.G.active_routing_server_url + "?z=" + OSRM.G.map.getZoom() + "&output=json&jsonp=%jsonp&alt=false&instructions=false";
	for(var i=0; i<chains[cur].length; ++i) {
		OSRM.G.markers.route[chains[cur][i]].marker.setOpacity(1);
		b += "&loc=" + 
		        OSRM.G.markers.route[chains[cur][i]].position.lat.toFixed(OSRM.C.PRECISION) + 
				"," + 
				OSRM.G.markers.route[chains[cur][i]].position.lng.toFixed(OSRM.C.PRECISION);
	}
	OSRM.JSONP.call(b, OSRM.RoutingGeometry.show, function(){}, OSRM.DEFAULTS.JSONP_TIMEOUT, "chain_draw")
	for(var i=1; i<chains[cur].length-1; ++i)
		OSRM.G.markers.route[chains[cur][i]].marker.setIcon(OSRM.G.icons['marker-via'])
	OSRM.G.markers.route[chains[cur][0]].marker.setIcon(OSRM.G.icons['marker-source'])
	OSRM.G.markers.route[chains[cur][chains[cur].length-1]].marker.setIcon(OSRM.G.icons['marker-target'])
	OSRM.G.markers.relabelViaMarkers()
	return false;
}
function showAllMarkers(){
	var b = OSRM.G.active_routing_server_url + "?z=" + OSRM.G.map.getZoom() + "&output=json&jsonp=%jsonp&alt=false&instructions=false";
	for(var i=0; i<glob.n; ++i) {
		OSRM.G.markers.route[i].marker.setOpacity(1);
		b += "&loc=" + 
		        OSRM.G.markers.route[i].position.lat.toFixed(OSRM.C.PRECISION) + "," + 
				OSRM.G.markers.route[i].position.lng.toFixed(OSRM.C.PRECISION);
	}
	OSRM.JSONP.call(b, OSRM.RoutingGeometry.show, function(){}, OSRM.DEFAULTS.JSONP_TIMEOUT, "chain_draw")
	for(var i=1; i<glob.n-1; ++i)
		OSRM.G.markers.route[i].marker.setIcon(OSRM.G.icons['marker-via'])
	OSRM.G.markers.route[0].marker.setIcon(OSRM.G.icons['marker-source'])
	OSRM.G.markers.route[glob.n-1].marker.setIcon(OSRM.G.icons['marker-target'])
	OSRM.G.markers.relabelViaMarkers()
	return false;
}
function showMarkersByPriority(l){
	var max=glob.counts[0][l];
	var min=glob.counts[0][l];
	for(var i=0; i<glob.n; ++i)
		if(min>glob.counts[i][l])
			min=glob.counts[i][l];
		else if(max<glob.counts[i][l])
			max=glob.counts[i][l];
	for(var i=0; i<glob.n; ++i)
	{
		OSRM.G.markers.route[i].marker.setOpacity((glob.counts[i][l]-min)/(max-min)*0.8+0.2);
		OSRM.G.markers.route[i].marker._icon.title=glob.counts[i][l];
	}
	return false;
}
// function showMarkersByPriorityDesc(){
	// OSRM.Routing.clearTransitMarkers()
	// priority=[];
	// var max=0;
	// for(var i=0; i<glob.n; ++i)
		// priority[i]=0;
	// for(var i=0;i<chains.length;i++)
		// for(var j=0; j<chains[i].length; ++j)
			// if(max<++priority[chains[i][j]])
				// max=priority[chains[i][j]];
	// var min=priority[0];
	// for(var i=0; i<glob.n; ++i)
		// if(min>priority[i])
			// min=priority[i];
	// for(var i=0; i<glob.n; ++i)
	// {
		// OSRM.G.markers.route[i].marker.setOpacity((1-(priority[i]-min)/(max-min))*0.8+0.2);
		// OSRM.G.markers.route[i].marker._icon.title=priority[i];
	// }
	// return false;
// }
function rgrand(min, max) {
    return Math.random() * (max - min) + min;
}
function randLatLng(p1, p2) {
    return L.latLng(rgrand(p1.lat, p2.lat), rgrand(p1.lng, p2.lng));
}
function addRandom(n) {
	for(var i=0;i<n; i++) {
		var index = OSRM.G.markers.setVia(0, randLatLng(OSRM.G.markers.route[0].position, 
			OSRM.G.markers.route[OSRM.G.markers.route.length-1].position));
		OSRM.G.markers.route[index].show();
	}
}

expand_attention = [];
nearest_expanded=[];
function expand_nearest_node(start, node)
{
	if(!expand_attention[node])
	{
		expand_attention[node]=true;
		nearest_expanded[start][node]=true;
		for(var i=0;i<glob.nearest[node].length;++i) 
			expand_nearest_node(start, glob.nearest[node][i])
	}
}
function expand_nearest() {
	nearest_expanded=[];
	nearest_bidir=[];
	for(var i=0;i<glob.n;++i) {	
		nearest_bidir[i]=[];
		nearest_expanded[i]=[];
		if(!nearest_expanded[0][i])
			nearest_expanded[0][i]=false;
		expand_attention = [];
		expand_nearest_node(i, i)
	}
	for(var i=0;i<glob.n;++i) {
		nearest_bidir[i].push(i)
		for(var j=0;j<i;++j)
			if(nearest_expanded[i][j] && nearest_expanded[j][i]) {
				nearest_bidir[i].push(j)
				nearest_bidir[j].push(i)
			}
	}
	glob.cores=[];
	for(var i=0;i<glob.n;++i) {
		nearest_bidir[i].sort(function(a,b){return a-b})
		b=true;
		for(var j=0;j<glob.cores.length;++j)
			if(glob.cores[j][0]==nearest_bidir[i][0]){
				b=false;
				break;
			}
		if(b) glob.cores.push(nearest_bidir[i])
	}
	glob.cores.sort(function(a,b){return b.length-a.length})
}
function showCores(){
	expand_nearest()
	var box=document.getElementById('clusters-box');
	var text='Выделить: <a href="#" onclick="return showAllMarkers()">все</a> <a href="#" onclick="return showMarkersByPriority()">значимые</a> <a href="#" onclick="return showMarkersByPriorityDesc()">редкие</a><ul>';
	for(var i=0;i<glob.cores.length;i++) {
		text+='<li><a href="#" onclick="return showCore('+i+')">'+i+'</a>=';
		text+=glob.cores[i].join('+');
		text+='</li>';
	}
	text+='</ul>';
	box.innerHTML=text;
}
function showCore(cur){
	for(var i=0; i<glob.n; ++i) {
		OSRM.G.markers.route[i].marker.setIcon(OSRM.G.icons['marker-highlight'])
		OSRM.G.markers.route[i].marker.setOpacity(0.8);
	}
	for(var i=0; i<glob.cores[cur].length; ++i)
		OSRM.G.markers.route[glob.cores[cur][i]].marker.setIcon(OSRM.G.icons['marker-via'])
	OSRM.G.markers.relabelViaMarkers()
	return false;
}
function compr(l){
	var sum=0;
	for(var i=0;i<glob.cores.length;i++) 
		if(glob.cores[i].length>=l)
			sum+=glob.cores[i].length
	return sum/glob.n*100;
}