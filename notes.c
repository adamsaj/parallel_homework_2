
particle particles[N];

struct box {
  particle* particles[10];
  int dx, dy;
};

box boxes[B];

//bin each point


for step in steps:

//part 2: apply forces
// For each box:
for(int i = 0; i < B; i++) {
  //For each point in the box:
  for(int p = 0; p < 10; p++) {
    //For every other point in the box:
    for(int p2 = 0; p2 < 10; p2++) {
      //Compute the force between the point p1 and the other point p2:
      apply_force( boxes[i].particles[p], boxes[i].particles[p2] );
    }
  }
  left =  box[i-1]; 
  right = box[i+1];
  top = box[i - box_size];
  bot = box[i + box_size];
  ll = box[i + box_size-1];

 }


for p in points:
  move(p);


for p in points:
  bin(p);

