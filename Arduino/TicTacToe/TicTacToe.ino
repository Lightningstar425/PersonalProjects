//Positive pins
int upper_left = 13;
int upper_mid = 12;
int upper_right = 11;

//Grounding pins
int left_up = 7;
int left_mid = 6;
int left_low = 5;

//Functions to light up each light

//Top left light
void upperLeftLight(){

  digitalWrite(upper_left, HIGH);
  digitalWrite(upper_mid, LOW);  
  digitalWrite(upper_right, LOW);
  
  digitalWrite(left_up, LOW);
  digitalWrite(left_mid, HIGH);
  digitalWrite(left_low, HIGH);
  
}

//Top mid Light
void upperMidLight(){
  digitalWrite(upper_left, LOW);
  digitalWrite(upper_mid, HIGH);  
  digitalWrite(upper_right, LOW);
  
  digitalWrite(left_up, LOW);
  digitalWrite(left_mid, HIGH);
  digitalWrite(left_low, HIGH);

}

//Top right light
void upperRightLight(){
  digitalWrite(upper_left, LOW);
  digitalWrite(upper_mid, LOW);  
  digitalWrite(upper_right, HIGH);
  
  digitalWrite(left_up, LOW);
  digitalWrite(left_mid, HIGH);
  digitalWrite(left_low, HIGH);
}

//Mid left light
void midLeftLight(){
  digitalWrite(upper_left, HIGH);
  digitalWrite(upper_mid, LOW);  
  digitalWrite(upper_right, LOW);
  
  digitalWrite(left_up, HIGH);
  digitalWrite(left_mid, LOW);
  digitalWrite(left_low, HIGH);
}

//middle light
void midMidLight(){
  digitalWrite(upper_left, LOW);
  digitalWrite(upper_mid, HIGH);  
  digitalWrite(upper_right, LOW);
  
  digitalWrite(left_up, HIGH);
  digitalWrite(left_mid, LOW);
  digitalWrite(left_low, HIGH);
}

//Middle right light
void midRightLight(){
  digitalWrite(upper_left, LOW);
  digitalWrite(upper_mid, LOW);  
  digitalWrite(upper_right, HIGH);
  
  digitalWrite(left_up, HIGH);
  digitalWrite(left_mid, LOW);
  digitalWrite(left_low, HIGH);
}

//Lower left light
void lowerLeftLight(){
  digitalWrite(upper_left, HIGH);
  digitalWrite(upper_mid, LOW);  
  digitalWrite(upper_right, LOW);
  
  digitalWrite(left_up, HIGH);
  digitalWrite(left_mid, HIGH);
  digitalWrite(left_low, LOW);
}

//Lower middle light
void lowerMidLight(){
  digitalWrite(upper_left, LOW);
  digitalWrite(upper_mid, HIGH);  
  digitalWrite(upper_right, LOW);
  
  digitalWrite(left_up, HIGH);
  digitalWrite(left_mid, HIGH);
  digitalWrite(left_low, LOW);
}

//Lower right light
void lowerRightLight(){
  digitalWrite(upper_left, LOW);
  digitalWrite(upper_mid, LOW);  
  digitalWrite(upper_right, HIGH);
  
  digitalWrite(left_up, HIGH);
  digitalWrite(left_mid, HIGH);
  digitalWrite(left_low, LOW);
}

void setup() {
  pinMode(upper_left,OUTPUT);
  pinMode(upper_mid,OUTPUT);
  pinMode(upper_right,OUTPUT);
  pinMode(left_up,OUTPUT);
  pinMode(left_mid,OUTPUT);
  pinMode(left_low,OUTPUT);

  digitalWrite(upper_left, LOW);
  digitalWrite(upper_mid, LOW);  
  digitalWrite(upper_right, LOW);
  
  digitalWrite(left_up, LOW);
  digitalWrite(left_mid, LOW);
  digitalWrite(left_low, LOW);
}

void loop() {
  
  upperLeftLight();
  delay(500);
  upperMidLight();
  delay(500);
  upperRightLight();
  delay(500);

  midLeftLight();
  delay(500);
  midMidLight();
  delay(500);
  midRightLight();
  delay(500);

  lowerLeftLight();
  delay(500);
  lowerMidLight();
  delay(500);
  lowerRightLight();
  delay(500);
  
}
