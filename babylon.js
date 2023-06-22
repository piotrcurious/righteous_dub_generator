// A function to generate a random bassline for 64 steps
function generateBassline() {
  // Initialize an array of 64 steps with random notes from C2 to G3
  let bassline = [];
  let notes = ["C2", "D2", "E2", "F2", "G2", "A2", "B2", "C3", "D3", "E3", "F3", "G3"];
  for (let i = 0; i < 64; i++) {
    let note = notes[Math.floor(Math.random() * notes.length)];
    bassline.push(note);
  }
  return bassline;
}

// A function to calculate the fitness of a bassline based on some criteria
function fitness(bassline) {
  // Initialize the fitness score to zero
  let score = 0;
  // Add some points for having a low note at the first step
  if (bassline[0] === "C2" || bassline[0] === "D2" || bassline[0] === "E2") {
    score += 10;
  }
  // Add some points for having a high note at the last step
  if (bassline[63] === "C3" || bassline[63] === "D3" || bassline[63] === "E3" || bassline[63] === "F3" || bassline[63] === "G3") {
    score += 10;
  }
  // Add some points for having a consistent rhythm
  let rhythm = true;
  for (let i = 0; i < 64; i += 4) {
    if (bassline[i] !== bassline[i + 1] || bassline[i + 1] !== bassline[i + 2] || bassline[i + 2] !== bassline[i + 3]) {
      rhythm = false;
      break;
    }
  }
  if (rhythm) {
    score += 20;
  }
  // Add some points for having some variation
  let variation = false;
  for (let i = 0; i < 64; i += 16) {
    if (bassline[i] !== bassline[i + 8] || bassline[i + 8] !== bassline[i + 16] || bassline[i + 16] !== bassline[i + 24]) {
      variation = true;
      break;
    }
  }
  if (variation) {
    score += 20;
  }
  
   // Add some points for having some harmony
   let harmony = true;
   for (let i = 0; i < 64; i +=4) {
     if (!isHarmonic(bassline[i],bassline[i+4])) {
       harmony = false;
       break;
     }
   }
   if (harmony) {
     score +=20;
   }

   // A helper function to check if two notes are harmonic
   function isHarmonic(note1, note2) {
     // A map of harmonic intervals
     let intervals = {
       "C2": ["C2", "G2", "C3"],
       "D2": ["D2", "A2", "D3"],
       "E2": ["E2", "B2", "E3"],
       "F2": ["F2", "C3", "F3"],
       "G2": ["G2", "D3", "G3"],
       "A2": ["A2", "E3"],
       "B2": ["B2", "F3"],
       "C3": ["C3", "G3"],
       "D3": ["D3"],
       "E3": ["E3"],
       "F3": ["F3"],
       "G3": ["G3"]
     };
     // Return true if note2 is in the interval list of note1
     return intervals[note1].includes(note2);
   }

   // Return the final fitness score
   return score;
}

// A function to create a new generation of basslines using crossover and mutation
function createGeneration(population, size) {
  // Initialize an array to store the new generation
  let generation = [];
  // Sort the population by fitness in descending order
  population.sort((a, b) => fitness(b) - fitness(a));
  // Select the top 10% of the population as the parents pool
  let pool = population.slice(0, Math.floor(size * 0.1));
  // Create new basslines by randomly picking two parents and crossing over their genes
  for (let i = 0; i < size; i++) {
    // Pick two random parents from the pool
    let parent1 = pool[Math.floor(Math.random() * pool.length)];
    let parent2 = pool[Math.floor(Math.random() * pool.length)];
    // Initialize an array to store the child bassline
    let child = [];
    // Loop through the 64 steps and randomly choose a gene from either parent
    for (let j = 0; j < 64; j++) {
      let gene = Math.random() < 0.5 ? parent1[j] : parent2[j];
      child.push(gene);
    }
    // Mutate the child bassline with a 2% chance of changing a random gene
    for (let j = 0; j < 64; j++) {
      if (Math.random() < 0.02) {
        let notes = ["C2", "D2", "E2", "F2", "G2", "A2", "B2", "C3", "D3", "E3", "F3", "G3"];
        let note = notes[Math.floor(Math.random() * notes.length)];
        child[j] = note;
      }
    }
    // Add the child bassline to the new generation
    generation.push(child);
  }
  // Return the new generation
  return generation;
}

// A function to create a game of life grid based on the user input
function createGrid(bassline) {
  // Initialize a 64 x 64 grid with all cells dead
  let grid = [];
  for (let i = 0; i < 64; i++) {
    let row = [];
    for (let j = 0; j < 64; j++) {
      row.push(0);
    }
    grid.push(row);
  }
  
   // Map each note in the bassline to a row index
   let noteToRow = {
     "C2": 0,
     "D2": 8,
     "E2": 16,
     "F2": 24,
     "G2": 32,
     "A2": 40,
     "B2": 48,
     "C3": 56,
     "D3": 57,
     "E3": 58,
     "F3": 59,
     "G3":60
   };

   // Loop through the bassline and make the corresponding cells alive
   for (let i =0; i<64;i++) {
     let note = bassline[i];
     let row = noteToRow[note];
     grid[row][i] =1;
   }

   // Return the grid
   return grid;
}

// A function to update the game of life grid based on the rules
function updateGrid(grid) {
  
   // Initialize a new grid with all cells dead
   let newGrid = [];
   for (let i=0;i<64;i++) {
     let row=[];
     for (let j=0;j<64;j++) {
       row.push(0);
     }
     newGrid.push(row);
   }

   // Loop through the old grid and apply the rules to each cell
   for (let i=0;i<64;i++) {
     for (let j=0;j<64;j++) {
       // Get the current state of the cell
       let state = grid[i][j];
       // Count the number of living neighbours
       let neighbours=0;
       for (let x=-1;x<=1;x++) {
         for (let y=-1;y<=1;y++) {
           if (x===0 && y===0) continue;
           if (i+x<0 || i+x>63 || j+y<0 || j+y>63) continue;
           neighbours +=grid[i+x][j+y];
         }
       }
       // Apply the rules based on the state and neighbours
       if (state ===1) { // living cell
         if (neighbours <2 || neighbours >3) { // underpopulation or overpopulation
           newGrid[i][j]=0; // die
         } else { // suitable environment
           newGrid[i][j]=1; // live
         }
       } else { // dead cell
         if (neighbours ===3) { // reproduction
           newGrid[i][j]=1; // live

         }
       }
     }
   }

   // Return the new grid
   return newGrid;
}

// A function to display the game of life grid on a canvas
function displayGrid(grid, canvas) {
  // Get the canvas context
  let ctx = canvas.getContext("2d");
  // Clear the canvas
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  // Set the cell size
  let cellSize = 10;
  // Loop through the grid and draw each cell
  for (let i = 0; i < 64; i++) {
    for (let j = 0; j < 64; j++) {
      // Get the state of the cell
      let state = grid[i][j];
      // Set the fill color based on the state
      if (state === 1) {
        ctx.fillStyle = "black";
      } else {
        ctx.fillStyle = "white";
      }
      // Draw a rectangle for the cell
      ctx.fillRect(j * cellSize, i * cellSize, cellSize, cellSize);
    }
  }
}

// A function to play the bassline on a web audio synth
function playBassline(bassline) {
  // Create a web audio context
  let audioCtx = new AudioContext();
  // Create an oscillator node for the synth
  let osc = audioCtx.createOscillator();
  // Set the oscillator type to square wave
  osc.type = "square";
  // Create a gain node for the volume control
  let gain = audioCtx.createGain();
  // Connect the oscillator to the gain node
  osc.connect(gain);
  // Connect the gain node to the destination (speakers)
  gain.connect(audioCtx.destination);
  // Set the initial volume to zero
  gain.gain.value = 0;
  // Start the oscillator
  osc.start();
  
   // Map each note in the bassline to a frequency in Hz
   let noteToFreq = {
     "C2":65.41,
     "D2":73.42,
     "E2":82.41,
     "F2":87.31,
     "G2":98.00,
     "A2":110.00,
     "B2":123.47,
     "C3":130.81,
     "D3":146.83,
     "E3":164.81,
     "F3":174.61,
     "G3":196.00
   };

   // Loop through the bassline and schedule each note to play at a certain time
   let stepDuration =0.25; // seconds per step
   for (let i=0;i<64;i++) {
     let note = bassline[i];
     let freq = noteToFreq[note];
     let time = audioCtx.currentTime + i * stepDuration;
     // Set the oscillator frequency to the note frequency at the start of the step
     osc.frequency.setValueAtTime(freq,time);
     // Ramp up the volume to 0.5 at the start of the step
     gain.gain.linearRampToValueAtTime(0.5,time);
     // Ramp down the volume to zero at the end of the step
     gain.gain.linearRampToValueAtTime(0,time + stepDuration);
   }

   // Stop the oscillator after playing all notes
   osc.stop(time + stepDuration);
}

// A function to run an evolutionary algorithm to optimize the bassline based on user input
function runEvolution(bassline, canvas) {
  
   // Initialize a population of basslines with random variations of the original bassline
   let population = [];
   let size =100; // population size
   for (let i=0;i<size;i++) {
     // Copy the original bassline
     let variation = bassline.slice();
     // Mutate some notes randomly with a 10% chance
     for (let j=0;j<64;j++) {
       if (Math.random()<0.1) {
         let notes = ["C2", "D2", "E2", "F2", "G2", "A2", "B2", "C3", "D3", "E3", "F3", "G3"];
         let note = notes[Math.floor(Math.random() * notes.length)];
         variation[j] = note;
       }
     }
     // Add the variation to the population
     population.push(variation);
   }

   // Initialize a generation counter
   let generation = 0;
   // Initialize a variable to store the best bassline
   let best = null;
   // Initialize a variable to store the best fitness
   let bestFitness = 0;
   // Initialize a variable to store the user input
   let userInput = null;

   // Create a function to update the user input based on mouse events on the canvas
   function updateUserInput(event) {
     // Get the mouse position relative to the canvas
     let rect = canvas.getBoundingClientRect();
     let x = event.clientX - rect.left;
     let y = event.clientY - rect.top;
     // Get the cell size
     let cellSize = 10;
     // Get the cell index based on the mouse position
     let i = Math.floor(y / cellSize);
     let j = Math.floor(x / cellSize);
     // Check if the mouse is inside the canvas and the cell is alive
     if (i >= 0 && i < 64 && j >= 0 && j < 64 && grid[i][j] === 1) {
       // Check if the user input is null or not
       if (userInput === null) {
         // Initialize the user input as an object with the cell index and the event type
         userInput = {
           i: i,
           j: j,
           type: event.type
         };
       } else {
         // Check if the user input has changed
         if (userInput.i !== i || userInput.j !== j || userInput.type !== event.type) {
           // Update the user input with the new cell index and event type
           userInput.i = i;
           userInput.j = j;
           userInput.type = event.type;
         }
       }
     }
   }

   // Add event listeners for mouse down, mouse up and mouse move on the canvas
   canvas.addEventListener("mousedown", updateUserInput);
   canvas.addEventListener("mouseup", updateUserInput);
   canvas.addEventListener("mousemove", updateUserInput);

   // Create a function to update the evolution based on user input and game of life rules
   function updateEvolution() {
     // Check if the user input is not null
     if (userInput !== null) {
       // Get the cell index and event type from the user input
       let i = userInput.i;
       let j = userInput.j;
       let type = userInput.type;
       // Get the note corresponding to the cell index
       let note = bassline[i];
       // Check if the event type is mousedown or mouseup
       if (type === "mousedown" || type === "mouseup") {
         // Loop through the population and modify each bassline based on the event type
         for (let k = 0; k < size; k++) {
           let variation = population[k];
           // If mousedown, increase the frequency of the note in the bassline by one step
           if (type === "mousedown") {
             variation[j] = increaseNote(note);
           }
           // If mouseup, decrease the frequency of the note in the bassline by one step
           if (type === "mouseup") {
             variation[j] = decreaseNote(note);
           }
         }
       }
       // Reset the user input to null
       userInput = null;
     }

     // Update the game of life grid based on the rules
     grid = updateGrid(grid);
     // Display the grid on the canvas
     displayGrid(grid, canvas);
     
      // A helper function to increase a note by one step in a chromatic scale
      function increaseNote(note) {
        // A map of notes and their next higher notes
        let nextNote = {
          "C2": "C#2",
          "C#2": "D2",
          "D2": "D#2",
          "D#2": "E2",
          "E2": "F2",
          "F2": "F#2",
          "F#2": "G2",
          "G2": "G#2",
          "G#2": "A2",
          "A2": "A#2",
          "A#2": "B2",
          "B2": "C3",
          "C3": "C#3",
          "C#3": "D3",
          "D3": "D#3",
          "D#3": "E3",
          "E3": "F3",
          "F3": "F#3",
          "F#3": "G3",
          "G3":"G#3"
        };
        // Return the next higher note or G#3 if already at highest note
        return nextNote[note] ||"G#3";
      }

      // A helper function to decrease a note by one step in a chromatic scale
      function decreaseNote(note) {
        // A map of notes and their next lower notes
        let prevNote = {
          "C2": "B1",
          "C#2": "C2",
          "D2": "C#2",
          "D#2": "D2",
          "E2": "D#2",
          "F2": "E2",
          "F#2": "F2",
          "G2": "F#2",
          "G#2": "G2",
          "A2": "G#2",
          "A#2": "A2",
          "B2": "A#2",
          "C3": "B2",
          "C#3": "C3",
          "D3": "C#3",
          "D#3": "D3",
          "E3": "D#3",
          "F3": "E3",
          "F#3": "F3",
          "G3":"F#3"
        };
        // Return the next lower note or B1 if already at lowest note
        return prevNote[note] ||"B1";
      }
     // Create a new generation of basslines using crossover and mutation
     population = createGeneration(population, size);
     // Increment the generation counter
     generation++;
     // Find the best bassline and its fitness in the current population
     let currentBest = null;
     let currentBestFitness = 0;
     for (let i = 0; i < size; i++) {
       let variation = population[i];
       let variationFitness = fitness(variation);
       if (variationFitness > currentBestFitness) {
         currentBest = variation;
         currentBestFitness = variationFitness;
       }
     }
     // Check if the best bassline has improved from the previous generation
     if (currentBestFitness > bestFitness) {
       // Update the best bassline and its fitness
       best = currentBest;
       bestFitness = currentBestFitness;
       // Play the best bassline
       playBassline(best);
     }
     // Log the generation number and the best fitness
     console.log("Generation: ", generation, ", Best fitness: ", bestFitness);
   }

   // Set an interval to update the evolution every 250 milliseconds
   let interval = setInterval(updateEvolution, 250);
}

// A function to make the game of life grid more interactive by allowing the user to toggle the state of a cell by clicking on it
function toggleCell(event) {
  // Get the mouse position relative to the canvas
  let rect = canvas.getBoundingClientRect();
  let x = event.clientX - rect.left;
  let y = event.clientY - rect.top;
  // Get the cell size
  let cellSize = 10;
  // Get the cell index based on the mouse position
  let i = Math.floor(y / cellSize);
  let j = Math.floor(x / cellSize);
  // Check if the mouse is inside the canvas
  if (i >= 0 && i < 64 && j >= 0 && j < 64) {
    // Toggle the state of the cell
    grid[i][j] = grid[i][j] === 1 ? 0 : 1;
    // Display the grid on the canvas
    displayGrid(grid, canvas);
    // Update the bassline based on the grid
    bassline = createBassline(grid);
    // Play the bassline
    playBassline(bassline);
  }
}

// Add an event listener for click on the canvas
canvas.addEventListener("click", toggleCell);

// A function to use a different crossover method that creates one offspring from many parents
function createGeneration(population, size) {
  // Initialize an array to store the new generation
  let generation = [];
  // Sort the population by fitness in descending order
  population.sort((a, b) => fitness(b) - fitness(a));
  // Select the top 10% of the population as the parents pool
  let pool = population.slice(0, Math.floor(size * 0.1));
  // Create new basslines by randomly picking many parents and crossing over their genes
  for (let i = 0; i < size; i++) {
    // Pick a random number of parents from the pool between 2 and 10
    let numParents = Math.floor(Math.random() * (10 - 2 + 1)) + 2;
    let parents = [];
    for (let j = 0; j < numParents; j++) {
      let parent = pool[Math.floor(Math.random() * pool.length)];
      parents.push(parent);
    }
    // Initialize an array to store the child bassline
    let child = [];
    // Loop through the 64 steps and randomly choose a gene from one of the parents
    for (let j = 0; j < 64; j++) {
      let gene = parents[Math.floor(Math.random() * numParents)][j];
      child.push(gene);
    }
    // Mutate the child bassline with a 2% chance of changing a random gene
    for (let j = 0; j < 64; j++) {
      if (Math.random() < 0.02) {
        let notes = ["C2", "D2", "E2", "F2", "G2", "A2", "B2", "C3", "D3", "E3", "F3", "G3"];
        let note = notes[Math.floor(Math.random() * notes.length)];
        child[j] = note;
      }
    }
    // Add the child bassline to the new generation
    generation.push(child);
  }
  // Return the new generation
  return generation;
}

// A function to make the game of life grid more colorful by using different colors for different notes
function displayGrid(grid, canvas) {
  // Get the canvas context
  let ctx = canvas.getContext("2d");
  // Clear the canvas
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  // Set the cell size
  let cellSize = 10;
  // Loop through the grid and draw each cell
  for (let i = 0; i < 64; i++) {
    for (let j = 0; j < 64; j++) {
      // Get the state of the cell
      let state = grid[i][j];
      // Set the fill color based on the state and the note
      if (state === 1) {
        let note = bassline[i];
        ctx.fillStyle = noteToColor(note);
      } else {
        ctx.fillStyle = "white";
      }
      // Draw a rectangle for the cell
      ctx.fillRect(j * cellSize, i * cellSize, cellSize, cellSize);
    }
  }

   // A helper function to map a note to a color
   function noteToColor(note) {
     // A map of notes and their colors
     let colors = {
       "C2": "red",
       "D2": "orange",
       "E2": "yellow",
       "F2": "green",
       "G2": "blue",
       "A2": "indigo",
       "B2": "violet",
       "C3": "pink",
       "D3": "brown",
       "E3": "gray",
       "F3": "black",
       "G3":"white"
     };
     // Return the color corresponding to the note or white if not found
     return colors[note] ||"white";
   }
}

// A function to make the game of life grid more dynamic by adding some randomness to the rules
function updateGrid(grid) {
  
   // Initialize a new grid with all cells dead
   let newGrid = [];
   for (let i=0;i<64;i++) {
     let row=[];
     for (let j=0;j<64;j++) {
       row.push(0);
     }
     newGrid.push(row);
   }

   // Loop through the old grid and apply the rules to each cell
   for (let i=0;i<64;i++) {
     for (let j=0;j<64;j++) {
       // Get the current state of the cell
       let state = grid[i][j];
       // Count the number of living neighbours
       let neighbours=0;
       for (let x=-1;x<=1;x++) {
         for (let y=-1;y<=1;y++) {
           if (x===0 && y===0) continue;
           if (i+x<0 || i+x>63 || j+y<0 || j+y>63) continue;
           neighbours +=grid[i+x][j+y];
         }
       }
       // Apply the rules based on the state and neighbours with some randomness
       if (state ===1) { // living cell
         if (neighbours <2 || neighbours >3) { // underpopulation or overpopulation
           newGrid[i][j]=0; // die with a 90% chance
           if (Math.random()<0.1) { // survive with a 10% chance
             newGrid[i][j]=1;
           }
         } else { // suitable environment
           newGrid[i][j]=1; // live with a 90% chance
           if (Math.random()<0.1) { // die with a 10% chance
             newGrid[i][j]=0;
           }
         }
       } else { // dead cell
         if (neighbours ===3) { // reproduction
           newGrid[i][j]=1; // live with a 90% chance
           if (Math.random()<0.1) { // stay dead with a 10% chance
             newGrid[i][j]=0;
           }
         } else { // no change
           newGrid[i][j]=0; // stay dead with a 90% chance
           if (Math.random()<0.1) { // live with a 10% chance
             newGrid[i][j]=1;

           }
         }
       }
     }
   }

   // Return the new grid
   return newGrid;
}

// A function to make the game of life grid more responsive by adjusting the cell size based on the canvas size
function displayGrid(grid, canvas) {
  // Get the canvas context
  let ctx = canvas.getContext("2d");
  // Clear the canvas
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  // Set the cell size based on the canvas size and the grid size
  let cellSize = Math.min(canvas.width, canvas.height) / 64;
  // Loop through the grid and draw each cell
  for (let i = 0; i < 64; i++) {
    for (let j = 0; j < 64; j++) {
      // Get the state of the cell
      let state = grid[i][j];
      // Set the fill color based on the state and the note
      if (state === 1) {
        let note = bassline[i];
        ctx.fillStyle = noteToColor(note);
      } else {
        ctx.fillStyle = "white";
      }
      // Draw a rectangle for the cell
      ctx.fillRect(j * cellSize, i * cellSize, cellSize, cellSize);
    }
  }

   // A helper function to map a note to a color
   function noteToColor(note) {
     // A map of notes and their colors
     let colors = {
       "C2": "red",
       "D2": "orange",
       "E2": "yellow",
       "F2": "green",
       "G2": "blue",
       "A2": "indigo",
       "B2": "violet",
       "C3": "pink",
       "D3": "brown",
       "E3": "gray",
       "F3": "black",
       "G3":"white"
     };
     // Return the color corresponding to the note or white if not found
     return colors[note] ||"white";
   }
}

// Add an event listener for window resize to adjust the canvas size
window.addEventListener("resize", function() {
  // Get the window width and height
  let width = window.innerWidth;
  let height = window.innerHeight;
  // Set the canvas width and height to be 80% of the window size
  canvas.width = width * 0.8;
  canvas.height = height * 0.8;
  // Display the grid on the canvas
  displayGrid(grid, canvas);
});

// A function to make the game of life grid more musical by changing the bassline based on the grid pattern
function createBassline(grid) {
  // Initialize an array of 64 steps with empty notes
  let bassline = [];
  for (let i = 0; i < 64; i++) {
    bassline.push("");
  }
  
   // Map each row index to a note
   let rowToNote = {
     0: "C2",
     8: "D2",
     16: "E2",
     24: "F2",
     32: "G2",
     40: "A2",
     48: "B2",
     56: "C3",
     57: "D3",
     58: "E3",
     59: "F3",
     60: "G3"
   };

   // Loop through the grid and find the first living cell in each column
   for (let j = 0; j < 64; j++) {
     for (let i = 0; i < 64; i++) {
       // Get the state of the cell
       let state = grid[i][j];
       // Check if the cell is alive
       if (state === 1) {
         // Get the note corresponding to the row index
         let note = rowToNote[i];
         // Set the bassline step to the note
         bassline[j] = note;
         // Break the loop
         break;
       }
     }
   }

   // Return the bassline
   return bassline;
}

// A function to play the bassline on a web audio synth with a different sound
function playBassline(bassline) {
  // Create a web audio context
  let audioCtx = new AudioContext();
  // Create an oscillator node for the synth
  let osc = audioCtx.createOscillator();
  // Set the oscillator type to sawtooth wave
  osc.type = "sawtooth";

  // Create a gain node for the volume control
  let gain = audioCtx.createGain();
  // Connect the oscillator to the gain node
  osc.connect(gain);
  // Connect the gain node to the destination (speakers)
  gain.connect(audioCtx.destination);
  // Set the initial volume to zero
  gain.gain.value = 0;
  // Start the oscillator
  osc.start();
  
   // Map each note in the bassline to a frequency in Hz
   let noteToFreq = {
     "C2":65.41,
     "D2":73.42,
     "E2":82.41,
     "F2":87.31,
     "G2":98.00,
     "A2":110.00,
     "B2":123.47,
     "C3":130.81,
     "D3":146.83,
     "E3":164.81,
     "F3":174.61,
     "G3":196.00
   };

   // Loop through the bassline and schedule each note to play at a certain time
   let stepDuration =0.25; // seconds per step
   for (let i=0;i<64;i++) {
     let note = bassline[i];
     let freq = noteToFreq[note];
     let time = audioCtx.currentTime + i * stepDuration;
     // Set the oscillator frequency to the note frequency at the start of the step
     osc.frequency.setValueAtTime(freq,time);
     // Ramp up the volume to 0.5 at the start of the step
     gain.gain.linearRampToValueAtTime(0.5,time);
     // Ramp down the volume to zero at the end of the step
     gain.gain.linearRampToValueAtTime(0,time + stepDuration);
   }

   // Stop the oscillator after playing all notes
   osc.stop(time + stepDuration);
}

// A function to make the game of life grid more complex by adding some predefined patterns
function createGrid(bassline) {
  // Initialize a 64 x 64 grid with all cells dead
  let grid = [];
  for (let i = 0; i < 64; i++) {
    let row = [];
    for (let j = 0; j < 64; j++) {
      row.push(0);
    }
    grid.push(row);
  }
  
   // Map each note in the bassline to a row index
   let noteToRow = {
     "C2": 0,
     "D2": 8,
     "E2": 16,
     "F2": 24,
     "G2": 32,
     "A2": 40,
     "B2": 48,
     "C3": 56,
     "D3": 57,
     "E3": 58,
     "F3": 59,
     "G3":60
   };

   // Loop through the bassline and make the corresponding cells alive
   for (let i =0; i<64;i++) {
     let note = bassline[i];
     let row = noteToRow[note];
     grid[row][i] =1;
   }

   // Add some predefined patterns to the grid based on the first note of the bassline
   let firstNote = bassline[0];
   switch (firstNote) {
     case "C2":
       // Add a glider pattern to the top left corner
       grid[1][0] = 1;
       grid[2][1] = 1;
       grid[0][2] = 1;
       grid[1][2] = 1;
       grid[2][2] = 1;
       break;
     case "D2":
       // Add a blinker pattern to the top left corner
       grid[1][0] = 1;
       grid[1][1] = 1;
       grid[1][2] = 1;
       break;
     case "E2":
       // Add a toad pattern to the top left corner
       grid[1][1] = 1;
       grid[1][2] = 1;
       grid[1][3] = 1;
       grid[2][0] = 1;
       grid[2][1] = 1;
       grid[2][2] = 1;
       break;
     case "F2":
       // Add a beacon pattern to the top left corner
       grid[0][0] = 1;
       grid[0][1] = 1;
       grid[1][0] = 1;
       grid[3][3] = 1;
       grid[3][4] = 1;
       grid[4][4] = 1;
       break;
     case "G2":
       // Add a pulsar pattern to the top left corner
       grid[0][2] = 1;
       grid[0][3] = 1;
       grid[0][4] = 1;
       grid[0][8] = 1;
       grid[0][9] = 1;
       grid[0][10] = 1;
       grid[2][0] = 1;
       grid[2][5] = 1;
       grid[2][7] = 1;
       grid[2][12] = 1;
       grid[3][2] = 1;
       grid[3][4] = 1;
       grid[3][8] = 1;
       grid[3][10] = 1;
       grid[4][2] = 1;
       grid[4][4] = 1;
       grid[4][8] = 1;
       grid[4][10] = 1;
       grid[5][0] = 1;
       grid[5][5] = 1;
       grid[5][7] = 1;
       grid[5][12] = 1;
       grid[7][0] = 1;
       grid[7][5] = 1;
       grid[7][7] = 1;
       grid[7][12] = 1;
       grid[8][2] = 1;
       grid[8][4] = 1;
       grid[8][8] = 1;
       grid[8][10] = 1;
       grid[9][2] = 1;
       grid[9][4] = 1;
       grid[9][8] = 1;
       grid[9][10] = 1;
       grid[10][0] = 1;
       grid[10][5] = 1;
       grid[10][7] = 1;
       grid[10][12] = 1;
       grid[12][2] = 1;
       grid[12][3] = 1;
       grid[12][4] = 1;
       grid[12][8] = 1;
       grid[12][9] = 1;
       grid[12][10] = 1; 
       break;  
     case "A2":
       // Add a pentadecathlon pattern to the top left corner
       grid[0][5] = 1; 
       grid[0][6] = 1; 
       grid[0][7] = 1; 
       grid[0][8] = 1; 
       grid[0][9] = 1; 
       grid[0][10] = 1; 
       grid[2][4] = 1; 
       grid[2][11] = 1; 
       grid[3][4] = 1; 
       grid[3][11] = 1; 
       grid[4][5] = 1; 
       grid[4][10] = 1; 
       grid[5][4] = 1; 
       grid[5][11] = 1; 
       grid[6][4] = 1; 
       grid[6][11] = 1; 
       grid[8][5] = 1; 
       grid[8][6] = 1; 
       grid[8][7] = 1; 
       grid[8][8] = 1; 
       grid[8][9] = 1; 
       grid[8][10] = 1; 
       break;
     case "B2":
       // Add a spaceship pattern to the top left corner
       grid[0][0] = 1;
       grid[0][3] = 1;
       grid[1][4] = 1;
       grid[2][0] = 1;
       grid[2][4] = 1;
       grid[3][1] = 1;
       grid[3][2] = 1;
       grid[3][3] = 1;
       grid[3][4] = 1;
       break;
     case "C3":
       // Add a glider gun pattern to the top left corner
       grid[0][24] = 1;
       grid[1][22] = 1;
       grid[1][24] = 1;
       grid[2][12] = 1;
       grid[2][13] = 1;
       grid[2][20] = 1;
       grid[2][21] = 1;
       grid[2][34] = 1;
       grid[2][35] = 1;
       grid[3][11] = 1;
       grid[3][15] = 1;
       grid[3][20] = 1;
       grid[3][21] = 1;
       grid[3][34] = 1;
       grid[3][35] = 1;
       grid[4][0] = 1;
       grid[4][1] = 1;
       grid[4][10] = 1;
       grid[4][16] = 1;
       grid[4][20] = 1;
       grid[4][21] = 1;
       grid[5][0] = 1;
       grid[5][1] = 1;
       grid[5][10] = 1;
       grid[5][14] = 1;
       grid[5][16] = 1;
       grid[5][17] = 1;
       grid[5][22] = 1;
       grid[5][24] = 1;
       grid[6][10] = 1;
       grid[6][16] = 1;
       grid[6][24] = 1;
       grid[7][11] = 1;
       grid[7][15] = 1;
       grid[8][12] = 1;
       grid[8][13] = 1; 
       break; 
     case "D3":
       // Add a r-pentomino pattern to the top left corner
       grid[0][0] = 1; 
       grid[0][1] = 1; 
       grid[0][2] = 1; 
       grid[1][0] = 1; 
       grid[2][1] = 1; 
       break; 
     case "E3":
       // Add a diehard pattern to the top left corner
       grid[0][6] = 1; 
       grid[1][0] = 1; 
       grid[1][1] = 1; 
       grid[2][1] = 1; 
       grid[2][5] = 1; 
       grid[2][6] = 1; 
       grid[2][7] = 1; 
       break; 
     case "F3":
       // Add an acorn pattern to the top left corner
       grid[0][0] = 1; 
       grid[0][2] = 1; 
       grid[2][3] = 1; 
       grid[2][4] = 1; 
       grid[2][5] = 1; 
       grid[2][6] = 1; 
       grid[4][1] = 1; 
       grid[4][2] = 1; 
       break; 
     case "G3":
       // Add a random pattern to the top left corner
       for (let i = 0; i < 8; i++) {
         for (let j = 0; j < 8; j++) {
           // Set the cell state to 1 with a 50% chance
           grid[i][j] = Math.random() < 0.5 ? 1 : 0;
         }
       }
       break;
   }

   // Return the grid
   return grid;
}
