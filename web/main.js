function getLatLon(address, callback) {
    console.log(`Loading coordinates for: ${address}`);
    let cachedData = localStorage.getItem(address);
    if (cachedData) {
        console.log(`Using cached data for: ${address}`);
        callback(JSON.parse(cachedData), true);
        return;
    }

    $.get("https://geocode.maps.co/search", { q: address, format: 'json' }, function(response) {
        if (response.length > 0) {
            let { lat, lon, type } = response[0];
            let data = { address, lat, lon, type, isCached: false };
            localStorage.setItem(address, JSON.stringify(data));
            console.log(`Coordinates loaded for: ${address}`);
            callback(data, false);
        } else {
            console.log(`No results found for: ${address}`);
            callback(null, false);
        }
    }).fail(function(xhr) {
        console.log(`Error loading coordinates for ${address}: ` + xhr.statusText);
        if (xhr.status === 429) setTimeout(() => getLatLon(address, callback), 1000);
    });
}

function calculateDistanceMatrix(coordinates, addresses, callback) {
    console.log("Generating distance matrix");
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
        console.log("Distance matrix generated");
        callback(distanceMatrix);
    }).fail(xhr => console.error("Error generating distance matrix: " + xhr.statusText));
}

function processAddresses(addresses, index, coordinates, callback) {
    console.log(`Processing address ${index + 1}/${addresses.length}`);
    if (index < addresses.length) {
        getLatLon(addresses[index], function(latlon, isCached) {
            if (latlon) coordinates.push({ lat: latlon.lat, lon: latlon.lon });
            let nextCall = () => processAddresses(addresses, index + 1, coordinates, callback);
            isCached ? nextCall() : setTimeout(nextCall, 1000);
        });
    } else {
        console.log("All addresses processed");
        callback(coordinates);
    }
}

Module.onRuntimeInitialized = _ => {
    let addresses = ["Queen's College, Taunton", "Musgrove Park Hospital", "Rumwell Farm Shop and Restaurant", "Taunton Train Station", "Museum of Somerset", "Creech St Michael Baptist Church", "Richard Huish College", "Trull Waterfall"];

    console.log("Starting address processing");
    processAddresses(addresses, 0, [], coordinates => {
        console.log("All coordinates:", coordinates);
        calculateDistanceMatrix(coordinates, addresses, distanceMatrix => {
            let csvContent = distanceMatrix.map(e => e.join(",")).join("\n");
            console.log("Starting TSP Solver");
            let maze = new Module.TSPSolver(csvContent, 100);
            console.log(csvContent);

            let routeString = maze.run();
            console.log("Route calculated:", routeString);

            // Process and URL encode the route string to create a URL for Google Maps
            let routeCities = routeString.split("\n")
                .map(line => line.trim().split(" -> ")[0])
                .map(city => encodeURIComponent(city));
            let url = "https://www.google.co.uk/maps/dir/" + routeCities.join('/');

            console.log("URL:", url);
        });
    });
};