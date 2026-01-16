// Clock face display pins
int pin_a = 3;
int pin_b = 5;
int pin_f = 4;
int pin_e = 11;
int pin_d = 2;
int pin_c = 9;
int pin_dp = 10;
int pin_g = 8;

//Change which thing is lighting up
int pin_1 = 13;
int pin_2 = 12;
int pin_3 = 6;
int pin_4 = 7;

void drawZero(){
  //Turn on
  digitalWrite(pin_a, HIGH);
  digitalWrite(pin_f, HIGH);
  digitalWrite(pin_b, HIGH);
  digitalWrite(pin_e, HIGH);
  digitalWrite(pin_c, HIGH);
  digitalWrite(pin_d, HIGH);

  //Turn off
  digitalWrite(pin_g, LOW);
}

void drawOne(){
  //Turn on
  digitalWrite(pin_b, HIGH);
  digitalWrite(pin_c, HIGH);

  //Turn off
  digitalWrite(pin_a, LOW);
  digitalWrite(pin_f, LOW);
  digitalWrite(pin_g, LOW);
  digitalWrite(pin_e, LOW);
  digitalWrite(pin_d, LOW);


}

void drawTwo(){
  //Turn on
  digitalWrite(pin_a, HIGH);
  digitalWrite(pin_b, HIGH);
  digitalWrite(pin_g, HIGH);
  digitalWrite(pin_e, HIGH);
  digitalWrite(pin_d, HIGH);

  //Turn off
  digitalWrite(pin_f, LOW);
  digitalWrite(pin_c, LOW);

}

void drawThree(){
  //Turn on
  digitalWrite(pin_a, HIGH);
  digitalWrite(pin_b, HIGH);
  digitalWrite(pin_g, HIGH);
  digitalWrite(pin_c, HIGH);
  digitalWrite(pin_d, HIGH);

  //Turn off
  digitalWrite(pin_f, LOW);
  digitalWrite(pin_e, LOW);


}

void drawFour(){
  //Turn on
  digitalWrite(pin_f, HIGH);
  digitalWrite(pin_g, HIGH);
  digitalWrite(pin_b, HIGH);
  digitalWrite(pin_c, HIGH);

  //Turn off
  digitalWrite(pin_a, LOW);
  digitalWrite(pin_e, LOW);
  digitalWrite(pin_d, LOW);

}

void drawFive(){
  //Turn on
  digitalWrite(pin_a, HIGH);
  digitalWrite(pin_f, HIGH);
  digitalWrite(pin_g, HIGH);
  digitalWrite(pin_c, HIGH);
  digitalWrite(pin_d, HIGH);

  //turn off
  digitalWrite(pin_b, LOW);
  digitalWrite(pin_e, LOW);
}

void drawSix(){
  digitalWrite(pin_a, HIGH);
  digitalWrite(pin_f, HIGH);
  digitalWrite(pin_g, HIGH);
  digitalWrite(pin_e, HIGH);
  digitalWrite(pin_c, HIGH);
  digitalWrite(pin_d, HIGH);

  digitalWrite(pin_b, LOW);

}

void drawSeven(){
  digitalWrite(pin_a, HIGH);
  digitalWrite(pin_b, HIGH);
  digitalWrite(pin_c, HIGH);

  digitalWrite(pin_f, LOW);
  digitalWrite(pin_g, LOW);
  digitalWrite(pin_e, LOW);
  digitalWrite(pin_d, LOW);

}

void drawEight(){
  digitalWrite(pin_a, HIGH);
  digitalWrite(pin_b, HIGH);
  digitalWrite(pin_c, HIGH);
  digitalWrite(pin_d, HIGH);
  digitalWrite(pin_e, HIGH);
  digitalWrite(pin_f, HIGH);
  digitalWrite(pin_g, HIGH);
}

void drawNine(){
  digitalWrite(pin_a, HIGH);
  digitalWrite(pin_b, HIGH);
  digitalWrite(pin_f, HIGH);
  digitalWrite(pin_g, HIGH);
  digitalWrite(pin_c, HIGH);

  digitalWrite(pin_e, LOW);
  digitalWrite(pin_d, LOW);
}

// List to call different functions
void (*drawFunctions[10])() = {
  drawZero, drawOne, drawTwo, drawThree, drawFour,
  drawFive, drawSix, drawSeven, drawEight, drawNine
};



void setup() {
  // put your setup code here, to run once:
  pinMode(pin_a, OUTPUT);
  pinMode(pin_b, OUTPUT);
  pinMode(pin_c, OUTPUT);
  pinMode(pin_d, OUTPUT);
  pinMode(pin_e, OUTPUT);
  pinMode(pin_f, OUTPUT);
  pinMode(pin_g, OUTPUT);
  pinMode(pin_dp, OUTPUT);

  pinMode(pin_1, OUTPUT);
  pinMode(pin_2, OUTPUT);
  pinMode(pin_3, OUTPUT);
  pinMode(pin_4, OUTPUT);


}

//ticker
int tick = 1000;
int counter = 0;

//Inital starts
int start_1 = 0;
int start_2 = 9;
int start_3 = 2;
int start_4 = 7;


void loop() {
  // put your main code here, to run repeatedly:  

  digitalWrite(pin_1, LOW);
  
  drawFunctions[(counter+start_1)%10]();
  delay(3);
  digitalWrite(pin_1, HIGH);

  digitalWrite(pin_2, LOW);
  drawFunctions[(counter+start_2)%10]();
  delay(3);
  digitalWrite(pin_2, HIGH);

  digitalWrite(pin_3, LOW);
  drawFunctions[(counter+start_3)%10]();
  delay(3);
  digitalWrite(pin_3, HIGH);


  digitalWrite(pin_4, LOW);
  drawFunctions[(counter+start_4)%10]();
  delay(3);
  digitalWrite(pin_4, HIGH);

  tick = tick + 1.0000000001;

  if (tick < 100){
    counter = 0;
  }
  else if (tick >= 100 && tick < 200){
      counter = 1;
  }
  else if (tick >= 200 && tick < 300){
      counter = 2;
  }
  else if (tick >= 300 && tick < 400){
      counter = 3;
  }
  else if (tick >= 400 && tick < 500){
      counter = 4;
  }
  else if (tick >= 500 && tick < 600){
      counter = 5;
  }
  else if (tick >= 600 && tick < 700){
      counter = 6;
  }
else if (tick >= 700 && tick < 800){
    counter = 7;
}
else if (tick >= 800 && tick < 900){
    counter = 8;
}
else {
    counter = 9;
}

if (tick > 1000){
  tick = 0;
}

  
}
