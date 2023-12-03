function getLatLon(address, callback) {
    let cachedData = localStorage.getItem(address);
    if (cachedData) {
        callback(JSON.parse(cachedData), true);
        return;
    }

    $.get("https://geocode.maps.co/search", { q: address, format: 'json' }, function(response) {
        if (response.length > 0) {
            let { lat, lon, type } = response[0];
            let data = { address, lat, lon, type, isCached: false };
            localStorage.setItem(address, JSON.stringify(data));
            callback(data, false);
        } else {
            console.log("No results found");
            callback(null, false);
        }
    }).fail(function(xhr) {
        console.log("Error: " + xhr.statusText);
        if (xhr.status === 429) setTimeout(() => getLatLon(address, callback), 1000);
    });
}

function calculateDistanceMatrix(coordinates, addresses, callback) {
    let matrixCacheKey = 'distanceMatrix_' + addresses.join('|');
    let cachedMatrix = localStorage.getItem(matrixCacheKey);

    if (cachedMatrix) {
        console.log("Using cached distance matrix");
        callback(JSON.parse(cachedMatrix));
        return;
    }

    let coordinateString = coordinates.map(coord => `${coord.lon},${coord.lat}`).join(';');
    $.get(`https://router.project-osrm.org/table/v1/driving/${coordinateString}`, { annotations: 'distance' }, function(response) {
        let distances = response.distances;
        let distanceMatrix = [[''].concat(addresses.map(address => `"${address}"`))];

        distances.forEach((row, index) => {
            let roundedRow = row.map(distance => Math.round(distance));
            distanceMatrix.push([`"${addresses[index]}"`].concat(roundedRow));
        });

        localStorage.setItem(matrixCacheKey, JSON.stringify(distanceMatrix));
        callback(distanceMatrix);
    }).fail(xhr => console.error("Error: " + xhr.statusText));
}

function processAddresses(addresses, index, coordinates, callback) {
    if (index < addresses.length) {
        getLatLon(addresses[index], function(latlon, isCached) {
            if (latlon) coordinates.push({ lat: latlon.lat, lon: latlon.lon });
            let nextCall = () => processAddresses(addresses, index + 1, coordinates, callback);
            isCached ? nextCall() : setTimeout(nextCall, 1000);
        });
    } else {
        callback(coordinates);
    }
}
Module.onRuntimeInitialized = _ => {
    // let addresses = ["Taunton", "St Ives", "Bristol", "London", "Middlesbrough", "Weston-super-Mare", "Birmingham", "Manchester", "Liverpool", "Newcastle", "Edinburgh", "Glasgow", "Cardiff", "Swansea", "Exeter", "Plymouth", "Bournemouth"];
    // let addresses = ["Los Angeles", "Kingman", "Grand Canyon Village", "Monument Valley", "Page", "Bryce Canyon City", "Zion National Park", "Las Vegas", "Death Valley", "Mammoth Lakes", "El Portal", "San Francisco", "Pismo Beach"];
    let addresses = ["Queen's College, Taunton", "Musgrove Park Hospital", "Rumwell Farm Shop and Restaurant", "Taunton Train Station", "Museum of Somerset", "Creech St Michael Baptist Church", "Richard Huish College", "Trull Waterfall"];


    processAddresses(addresses, 0, [], coordinates => {
        console.log("All coordinates:", coordinates);
        calculateDistanceMatrix(coordinates, addresses, distanceMatrix => {
            let csvContent = distanceMatrix.map(e => e.join(",")).join("\n");
            let maze = new Module.TSPSolver(csvContent, 100);
            console.log(csvContent);

            let routeString = maze.run();
            console.log("Route:", routeString);

            // Process and URL encode the route string to create a URL for Google Maps
            let routeCities = routeString.split("\n")
                .map(line => line.trim().split(" -> ")[0])
                .map(city => encodeURIComponent(city));
            let url = "https://www.google.co.uk/maps/dir/" + routeCities.join('/');

            console.log("URL:", url);
        });
    });
};