function getLatLon(address, callback) {
    console.log(`Loading coordinates for: ${address}`);
    let cachedData = localStorage.getItem(address);
    if (cachedData) {
        console.log(`Using cached data for: ${address}`);
        callback(JSON.parse(cachedData), true);
        return;
    }

    $.get("https://nominatim.openstreetmap.org/search.php", { q: address, format: 'json' }, function(response) {
        if (response.length > 0) {
            let { lat, lon, type } = response[0];
            let data = { address, lat, lon, type, isCached: false };
            localStorage.setItem(address, JSON.stringify(data));
            console.log(`Coordinates loaded for: ${address}`);
            callback(data, false);
        } else {
            alert(`No results found for: ${address}`);
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
    // Progress bar starts from 5% and ends at 50%
    UpdateStatus(`Processing address ${addresses[index]} (${index + 1}/${addresses.length})`, 5 + (45 * index / addresses.length));
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

var module = null;

function UpdateStatus(status, progress) {
    document.getElementById('status').innerHTML = status;
    document.getElementById('statusBar').style.width = progress + '%';
}

document.getElementById('demoRoutes').addEventListener('click', function(e) {
    let addresses = e.target.dataset.demoroute.split(',').map(city => city.trim());
    document.getElementById('cityInput').value = addresses.join('\n');
});

document.getElementById('resetCache').addEventListener('click', function(e) {
    // input for confirmation
    let input = prompt("Are you sure you want to clear the cache? This will cause all coordinates to be reloaded from the API. Type 'yes' to confirm.");
    if (input === "yes") {
        localStorage.clear();
        alert("Cache Cleared!");
    }
});

Module.onRuntimeInitialized = _ => {
    console.log("Module loaded");

    document.getElementById('inputContainer').addEventListener('submit', function(e) {
        e.preventDefault();

        let input = document.getElementById('cityInput').value;
        let addresses = input.split('\n').map(city => city.trim());
        addresses = addresses.filter((city, index, self) => city.length > 0 && self.indexOf(city) === index);
        UpdateStatus("Processing addresses...", 5);
        CalculateRoute(addresses);
    });

    function CalculateRoute(addresses) {
        processAddresses(addresses, 0, [], coordinates => {
            UpdateStatus("Loading distance matrix...", 60);
            calculateDistanceMatrix(coordinates, addresses, distanceMatrix => {
                let csvContent = distanceMatrix.map(e => e.join(",")).join("\n");
                UpdateStatus("Calculating route...", 80);
                let solver = new Module.TSPSolver(csvContent, 100);
                console.log(csvContent);
    
                solver.run();

                let routeString = solver.GetBestRoute();

                UpdateStatus("Finishing up...", 90);

                // Process and URL encode the route string to create a URL for Google Maps
                let routeCities = routeString.split("\n")
                    .map(line => line.trim().replace(/"/g, ''))
                    .map(city => encodeURIComponent(city));
                let url = "https://www.google.co.uk/maps/dir/" + routeCities.join('/');
        
                console.log("URL:", url);

                UpdateStatus("Complete!", 100);
        
                document.getElementById('routeDisplay').innerHTML = routeString.split("\n").map(line => line.trim().replace(/"/g, '')).filter(line => line).join(' ➡️ ');
                document.getElementById('distanceDisplay').innerHTML = "Distance: " + Math.round(solver.GetBestRouteLength() / 1609.344 * 10) / 10 + " miles";
                document.getElementById('routeLink').setAttribute('href', url);
                document.getElementById('routeLink').innerHTML = 'View on Google Maps';
            });
        });
    }
};
