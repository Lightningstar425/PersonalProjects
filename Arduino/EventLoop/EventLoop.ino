//Libraries
#include <PriorityQueue.h>
#include <QueueList_Modified.h>


//Define waitTime
long waitTime = 5000;
long action = 0;

// Event loop times
long currentTime = 0;

//Pin colors
int greenPin = 8;
long greenWait = 3000;

int redPin =  7;
long redWait = 5000;

typedef struct Event {
  long execute_at_time;
  void (*run)(Event *event);
};

//Loops and ques
Event *event;
PriorityQueue<Event*> *queue;

void hello(Event *event){
  Serial.println("hello");
  event = new Event();
  event->execute_at_time = waitTime + millis();
  event->run = hello;
  queue->push(event);

}

void greenLight(Event *event){
  digitalWrite(greenPin, !digitalRead(greenPin));
  Serial.println("green");
  event = new Event();
  event->execute_at_time = greenWait + millis();
  event->run = greenLight;
  queue->push(event);

}

void redLight(Event *event){
  digitalWrite(redPin, !digitalRead(redPin));
  Serial.println("red");
  event = new Event();
  event->execute_at_time = redWait + millis();
  event->run = redLight;
  queue->push(event);

}


bool compareTime(Event *a, Event *b){
  return a->execute_at_time < b->execute_at_time;
}

void setup() {
  Serial.begin(9600);
  queue = new PriorityQueue<Event*>(compareTime);


  event = new Event();
  event->execute_at_time = redWait;
  event->run = redLight;
  queue->push(event);

  event = new Event();
  event->execute_at_time = greenWait;
  event->run = greenLight;
  queue->push(event);
}

void loop() {
  currentTime = millis();
  if (!queue->isEmpty()){
    action = queue->peek()->execute_at_time;
    if (action <= currentTime){
      event = queue->pop();
      event->run(event);
      delete event;
    }
  }
}
