

#include <AccelStepper.h>

//Motor pins
int pin4 = 11;
int pin3 = 10;
int pin2 = 9;
int pin1 = 8;

int MotorInterfaceType = 8;

AccelStepper stepper = AccelStepper(MotorInterfaceType, pin1, pin3, pin2, pin4);
int spr = 4800;

void setup() {
  // put your setup code here, to run once:
  stepper.setMaxSpeed(1000);
  stepper.setAcceleration(200);
  stepper.setSpeed(600);
}

void loop() {
  // put your main code here, to run repeatedly:
  //stepper.moveTo(stepper.currentPosition() + spr);
  //stepper.runToPosition();
  stepper.runSpeed();
}
