//Libraries
#include <PriorityQueue.h>
#include <QueueList_Modified.h>
#include <AccelStepper.h>

//Event Structure
typedef struct Event {
  int id;
  unsigned long long timeToRun;
  void (*run)(Event *event);
};

//Initialize pointers
PriorityQueue<Event *> *queue;

//Ultrasonic sensor pins
const int trigPin = 3;
const int echoPin = 4;
int lightPin = 6;
uint8_t lightOn = LOW;
float duration, distance;


//Motor pins
int pin4 = 11;
int pin3 = 10;
int pin2 = 9;
int pin1 = 8;
int motorPower = 12;

int MotorInterfaceType = 8;
AccelStepper stepper = AccelStepper(MotorInterfaceType, pin1, pin3, pin2, pin4);

const long spr = 2048;
int count = 0;
float past_dis;


//Time delays (in Microsecounds)
const unsigned long firstTrigDelay = 2;
const unsigned long secoundTrigDelay = 15;
const unsigned long timeBetweenEchos = 1000000;

void lightSwitch(Event *event) {
  digitalWrite(lightPin, !digitalRead(lightPin));
}

void setMotor(Event *event) {
  stepper.moveTo(stepper.currentPosition() + 2048);
}

void firstTrigLow(Event *event) {
  Event *next = new Event();
  next->timeToRun = firstTrigDelay + micros();
  next->run = trigHigh;
  queue->push(next);
  digitalWrite(trigPin, LOW);
}

void trigHigh(Event *event) {
  Event *next = new Event();
  next->timeToRun = secoundTrigDelay + micros();
  next->run = lastTrigLow;
  queue->push(next);
  digitalWrite(trigPin, HIGH);
}

void lastTrigLow(Event *event) {

  digitalWrite(trigPin, LOW);
  //Serial.println(micros());
  duration = pulseIn(echoPin, HIGH);
  //Serial.println(micros());
  past_dis = distance;
  distance = (duration * 0.0343) / 2;

  Serial.println(distance);
  if (distance < 5 && past_dis < 5) {
    count++;
    if (count >= 5) {
      Event *motorEvent = new Event();
      motorEvent->timeToRun = micros();
      motorEvent->run = setMotor;
      queue->push(motorEvent);
      count = 0;
    }


  } else {
    count = 0;
  }

  Event *next = new Event();
  next->timeToRun = timeBetweenEchos + micros();
  next->run = firstTrigLow;
  queue->push(next);
}

bool compareTime(Event *a, Event *b) {
  return a->timeToRun < b->timeToRun;
}

void setup() {
  //Initalize pinss
  pinMode(trigPin, OUTPUT);
  pinMode(echoPin, INPUT);
  pinMode(lightPin, OUTPUT);
  pinMode(motorPower, OUTPUT);
  digitalWrite(motorPower, HIGH);
  Serial.begin(9600);
  stepper.setMaxSpeed(1000);
  stepper.setAcceleration(200);
  stepper.setSpeed(600);

  //Start Queue and echoing
  queue = new PriorityQueue<Event *>(compareTime);
  Event *event = new Event();
  event->timeToRun = 0;
  event->run = firstTrigLow;
  queue->push(event);
}

void loop() {

  unsigned long long currentTime = micros();
  if (!queue->isEmpty()) {
    unsigned long long actionTime = queue->peek()->timeToRun;
    if (actionTime <= currentTime) {
      Event *current = queue->pop();
      current->run(current);
      delete current;
    }
  }
  if (distance < 5 && past_dis > 5) {
    Event *event = new Event();
    event->timeToRun = currentTime;
    event->run = lightSwitch;
    queue->push(event);
  }


  stepper.run();
}